# Bad controls: mediators and colliders

A regression that "controls for" many variables is not better for it. A control helps when it blocks a back-door path,
costs at most precision when it lies on no path, and does harm when it lies on the causal path or when conditioning on it
opens a closed path. Angrist and Pischke (2009) call the last two kinds **bad controls**. This section computes the damage
in the linear models of [Figure 25.5.1](#fig-cau-bad-dags) and describes the related bias from selecting the sample on a
collider.

\begin{center}
\begin{tikzpicture}[>=Stealth, line cap=round]
\definecolor{ink}{HTML}{1B1F27}\definecolor{accent}{HTML}{2F5F96}\definecolor{draft}{HTML}{A45C25}
\tikzset{var/.style={draw=ink, circle, minimum size=7.5mm, inner sep=0pt, font=\small}}
\node[var] (X1) at (0,0) {$X$}; \node[var, fill=draft!15] (C1) at (1.3,1.0) {$C$}; \node[var] (Y1) at (2.6,0) {$Y$};
\draw[->, thick, accent] (X1) -- node[above left, font=\small] {$\delta$} (C1);
\draw[->, thick, accent] (C1) -- node[above right, font=\small] {$\gamma$} (Y1);
\draw[->, thick, accent] (X1) -- node[below, font=\small] {$\theta$} (Y1);
\node[font=\small] at (1.3,-0.8) {(a) mediator};
\node[var] (X2) at (3.9,0) {$X$}; \node[var, fill=draft!15] (C2) at (5.2,1.0) {$C$}; \node[var] (Y2) at (6.5,0) {$Y$};
\draw[->, thick, ink] (X2) -- node[above left, font=\small] {$a$} (C2);
\draw[->, thick, ink] (Y2) -- node[above right, font=\small] {$b$} (C2);
\draw[->, thick, accent] (X2) -- node[below, font=\small] {$\tau$} (Y2);
\node[font=\small] at (5.2,-0.8) {(b) common effect};
\node[var, dashed] (H1) at (7.8,1.6) {$H_1$}; \node[var, dashed] (H2) at (10.4,1.6) {$H_2$};
\node[var, fill=draft!15] (C3) at (9.1,0.9) {$C$};
\node[var] (X3) at (7.8,0) {$X$}; \node[var] (Y3) at (10.4,0) {$Y$};
\draw[->, thick, draft] (H1) -- node[left, font=\small] {$a$} (X3);
\draw[->, thick, draft] (H1) -- node[above, font=\small] {$c_1$} (C3);
\draw[->, thick, draft] (H2) -- node[above, font=\small] {$c_2$} (C3);
\draw[->, thick, draft] (H2) -- node[right, font=\small] {$b$} (Y3);
\draw[->, thick, accent] (X3) -- node[below, font=\small] {$\tau$} (Y3);
\node[font=\small] at (9.1,-0.8) {(c) M-structure};
\end{tikzpicture}
\end{center}

[**Figure 25.5.1.** Three bad controls \( C \) (shaded). (a) \( C \) is a mediator. (b) \( C \) is a common effect of \( X \) and \( Y \).
(c) \( C \) is measured before \( X \) but is a collider on the back-door path through the unrecorded \( H_1 \) and
\( H_2 \).]{#fig-cau-bad-dags}

## Three computations

::: {#thm-cau-bad-controls}
[Bad controls in linear models]

In each model below the disturbances \( U_X,U_Y,U_C \) (and, in (c), the unrecorded variables \( H_1,H_2 \)) are independent with mean zero and variances
\( \sigma_X^2,\sigma_Y^2,\sigma_C^2>0 \) (and \( 1 \)). Let \( \tau \) be the total effect of \( X \) on \( Y \). In every case the slope of
\( L(Y\mid X) \) equals \( \tau \), so the simple regression is consistent, while the coefficient \( \beta_{X\mid C} \) of \( X \) in
\( L(Y\mid X,C) \) is as follows.

::: {.enumerate options="label=(\alph*)"}
1. **Mediator:** \( X=U_X \), \( C=\delta X+U_C \), \( Y=\theta X+\gamma C+U_Y \). Here \( \tau=\theta+\gamma\delta \) and
   \( \beta_{X\mid C}=\theta \): adjustment removes the indirect effect \( \gamma\delta \).

2. **Common effect:** \( X=U_X \), \( Y=\tau X+U_Y \), \( C=aX+bY+U_C \). Then
   \[
   \beta_{X\mid C}=\tau-\frac{(a+b\tau)\,b\,\sigma_Y^2}{b^2\sigma_Y^2+\sigma_C^2}.
   \]{#eq-cau-bias-common-effect}

3. **M-structure:** \( X=aH_1+U_X \), \( C=c_1H_1+c_2H_2+U_C \), \( Y=\tau X+bH_2+U_Y \). Then
   \[
   \begin{aligned}
   \beta_{X\mid C}&=\tau-\frac{a\,b\,c_1c_2}{D},\\
   D&=(a^2+\sigma_X^2)(c_1^2+c_2^2+\sigma_C^2)-a^2c_1^2>0 .
   \end{aligned}
   \]{#eq-cau-bias-m}

:::

:::

::: {.proof}
(a) \( Y-\theta X-\gamma C=U_Y \) is uncorrelated with \( X \) and \( C \), so \( L(Y\mid X,C)=\theta X+\gamma C \). Without \( C \),
\( Y=(\theta+\gamma\delta)X+(\gamma U_C+U_Y) \), and the bracket is uncorrelated with \( X \), so the slope is \( \tau \).

(b) Without \( C \), \( U_Y \) is uncorrelated with \( X \) and the slope is \( \tau \). With \( C \), write
\( E=C-(a+b\tau)X=bU_Y+U_C \), which is uncorrelated with \( X \). The span of \( 1,X,C \) is the span of \( 1,X,E \), with \( E \)
orthogonal to \( 1 \) and \( X \), so \( L(U_Y\mid X,C)=L(U_Y\mid X)+L(U_Y\mid E)=0+\lambda E \), where
\( \lambda=\Cov(U_Y,E)/\Var(E)=b\sigma_Y^2/(b^2\sigma_Y^2+\sigma_C^2) \). Hence
\( L(Y\mid X,C)=\tau X+\lambda C-\lambda(a+b\tau)X \), whose coefficient of \( X \) is @eq-cau-bias-common-effect.

(c) Without \( C \), \( bH_2+U_Y \) is uncorrelated with \( X \) and the slope is \( \tau \). With \( C \),
\( L(Y\mid X,C)=\tau X+b\,L(H_2\mid X,C) \), since \( U_Y \) is uncorrelated with \( X \) and \( C \). The coefficients of
\( L(H_2\mid X,C) \) solve the population equations of @thm-proj-blp, with covariance matrix of \( (X,C) \) equal to
\[
\begin{pmatrix}a^2+\sigma_X^2&ac_1\\ac_1&c_1^2+c_2^2+\sigma_C^2\end{pmatrix}
\]
and right-hand side \( (\Cov(X,H_2),\Cov(C,H_2))=(0,c_2) \). The determinant is \( D \), positive since the matrix is a
positive definite covariance matrix, and the first coefficient is \( -ac_1c_2/D \). Multiply by \( b \).
:::

In (a) the adjusted coefficient is the **direct effect**, with the mediator held fixed by intervention. That may be of
interest, but only if nothing confounds the mediator and the outcome (@exr-cau-mediator-confounding); a policy decision
needs \( \tau \) ([Section 18.1](../ch18-covariance-and-design/01-ancova.html) gave the experimental version).

In (b) the bias exists even when \( \tau=0 \): adjusting for a consequence of both variables produces an "effect"
\( -ab\sigma_Y^2/(b^2\sigma_Y^2+\sigma_C^2) \) where there is none, negative when \( a \) and \( b \) have the same sign. This is the
collider of @prp-cau-three-structures(c).

In (c) the control is measured before \( X \), so the advice to adjust only for pre-treatment variables does not exclude it,
yet adjustment opens the back-door path \( X\leftarrow H_1\to C\leftarrow H_2\to Y \). This **M-bias**, named after the
graph's shape, was pointed out by Greenland, Pearl and Robins (1999). It is a product of four path coefficients over \( D \)
and tends to be small unless all four are large (Ding and Miratrix 2015). When the collider is also a confounder, with
arrows \( C\to X \) and \( C\to Y \), neither choice removes the bias; Ding and Miratrix call this butterfly bias.

::: {#exm-cau-bad-controls}
[Three bad controls by simulation]

Take, with unit variances throughout:

- a mediator with \( \theta=0.5 \), \( \delta=1.2 \), \( \gamma=1 \), so that \( \tau=1.7 \) and \( \beta_{X\mid C}=0.5 \);
- a common effect with \( \tau=1 \), \( a=0.5 \), \( b=1 \), for which @eq-cau-bias-common-effect gives
  \( \beta_{X\mid C}=1-(1.5)(1)/2=0.25 \);
- an M-structure with \( \tau=1 \) and \( a=b=c_1=c_2=1.5 \), for which
  \( D=\tfrac{13}{4}\cdot\tfrac{11}2-\tfrac{81}{16}=\tfrac{205}{16}=12.8125 \) and
  \( \beta_{X\mid C}=1-\tfrac{81}{205}=0.605 \).

Over \( 4000 \) samples of \( n=400 \), the least squares coefficients of \( X \) without and with \( C \) average
\( 1.700 \) and \( 0.499 \) for the mediator, \( 1.000 \) and
\( 0.251 \) for the common effect, and \( 0.999 \) and
\( 0.605 \) for the M-structure ([Figure 25.5.2](#fig-cau-bad-controls)). In the M-structure the
adjusted estimate is also slightly *less* variable (standard deviation \( 0.046 \) against
\( 0.050 \)), a reminder that a smaller standard error says nothing about bias.
:::

::: {when-format="html"}
![**Figure 25.5.2.** Sampling distributions of the least squares coefficient of \( X \) without (blue) and with (brown)
adjustment for \( C \), in the three models of @exm-cau-bad-controls; \( 4000 \) samples of size \( 400 \). The dashed line
is the total effect \( \tau \).](bad_controls.svg){#fig-cau-bad-controls width=100%}
:::

::: {when-format="pdf"}
![Sampling distributions of the least squares coefficient of \( X \) without (blue) and with (brown)
adjustment for \( C \), in the three models of @exm-cau-bad-controls; \( 4000 \) samples of size \( 400 \). The dashed line
is the total effect \( \tau \).](bad_controls.pdf){width=100%}
:::

```{.python .run #cell-bad-controls-models}
import numpy as np

rng = np.random.default_rng(2505)

def mediator(n):          # X -> C -> Y and X -> Y; total effect 0.5 + 1.2 * 1.0 = 1.7
    x = rng.normal(size=n); c = 1.2 * x + rng.normal(size=n)
    return x, c, 0.5 * x + 1.0 * c + rng.normal(size=n)

def collider(n):          # X -> Y, and X -> C <- Y; total effect 1
    x = rng.normal(size=n); y = 1.0 * x + rng.normal(size=n)
    return x, 0.5 * x + 1.0 * y + rng.normal(size=n), y

def m_bias(n):            # X <- H1 -> C <- H2 -> Y, and X -> Y; total effect 1
    h1, h2 = rng.normal(size=(2, n))
    x = 1.5 * h1 + rng.normal(size=n)
    c = 1.5 * h1 + 1.5 * h2 + rng.normal(size=n)
    return x, c, 1.0 * x + 1.5 * h2 + rng.normal(size=n)

def slopes(x, c, y):
    """Coefficient of x without and with adjustment for c."""
    one = np.ones(len(x))
    b_short = np.linalg.lstsq(np.column_stack([one, x]), y, rcond=None)[0][1]
    b_long = np.linalg.lstsq(np.column_stack([one, x, c]), y, rcond=None)[0][1]
    return b_short, b_long
```

```{.python .run #cell-bad-controls-quick}
n, reps = 400, 500
for name, model in [("mediator", mediator), ("collider", collider), ("M-bias", m_bias)]:
    est = np.array([slopes(*model(n)) for _ in range(reps)])
    print(f"{name:9s} mean without C {est[:, 0].mean():.3f}   mean with C {est[:, 1].mean():.3f}")
```

## Post-treatment variables in general

A variable affected by \( X \) removes part of the effect if it lies on a directed path to \( Y \), and opens a path if \( Y \)
or a cause of \( Y \) also affects it; its descendants act partially (@exr-cau-collider-descendant). A descendant of \( X \)
with no other connection to \( Y \) costs only precision (@exr-cau-descendant-only). "Adjust only for pre-treatment
variables" is sound in randomized experiments ([Section 25.6](06-randomized-adjustment.html)); in observational data it is
neither sufficient, by the M-structure, nor necessary, since a later measurement can proxy an unrecorded confounder. The
graph, not the timing, decides. Adjusting for a pure cause of \( X \) can even *increase* the bias from an unrecorded
confounder (@exr-cau-bias-amplification); Pearl (2010) calls such variables bias amplifiers.

## Selection on a collider

A sample selected on a common effect of \( X \) and \( Y \) is conditioned on a collider without any adjustment. Berkson
(1946) saw two diseases, independent in the population, associated among hospital patients, since either brings a patient to
hospital; Elwert and Winship (2014) survey this *endogenous selection bias*. In the common-effect model of
@exm-cau-bad-controls, observing only the units with \( C>1 \) (about \( 0.31 \) of them) moves the slope of \( Y \) on
\( X \) from \( 1.002 \) to \( 0.584 \).

```{.python .run #cell-bad-controls-selection}
x, c, y = collider(200_000)
keep = c > 1                                          # only units with a large C are observed
slope_all = np.polyfit(x, y, 1)[0]
slope_sel = np.polyfit(x[keep], y[keep], 1)[0]
print(f"slope of Y on X: all units {slope_all:.3f}, units with C > 1 {slope_sel:.3f}")
```

::: {.warning}
A variable belongs in a causal regression because of its position in the graph, not because it improves \( R^2 \) or is
"significant"; predictive model selection ([Chapter 29](../ch29-model-selection/index.html)) readily selects mediators and colliders.
:::

## Exercises

### A. Check your understanding

::: {#exr-cau-classify}
[A1]

A study of the effect of a training programme (\( X \)) on earnings two years later (\( Y \)) considers adjusting for: (i)
education before the programme; (ii) the job obtained six months after the programme; (iii) whether the person was still
living in the region at the end of the study, which depends on both training and earnings; (iv) the day of the week on which the
person applied. Classify each as a confounder, mediator, collider or irrelevant, under a plausible graph, and say which
should be adjusted for.
:::

### B. Practice

::: {#exr-cau-descendant-only}
[B1]

Let \( X=U_X \), \( Y=\tau X+U_Y \) and \( D=dX+U_D \), where \( D \) is affected by \( X \) but has no connection with \( Y \). Show that
the coefficient of \( X \) in \( L(Y\mid X,D) \) is \( \tau \), and that the asymptotic variance of the least squares coefficient
is larger with \( D \) than without it by the factor \( (d^2\sigma_X^2+\sigma_D^2)/\sigma_D^2 \).
:::

::: {#exr-cau-mediator-confounding}
[B2]

In the mediator model add an unrecorded \( U \), with variance \( \sigma_U^2 \), that affects both the mediator and the outcome:
\( C=\delta X+gU+U_C \), \( Y=\theta X+\gamma C+hU+U_Y \). Show that
\[
\beta_{X\mid C}=\theta-\frac{\delta\,g\,h\,\sigma_U^2}{g^2\sigma_U^2+\sigma_C^2},
\]
so that adjusting for the mediator does not even recover the direct effect \( \theta \).
:::

::: {.solution}
\( Y=\theta X+\gamma C+hU+U_Y \) and \( U_Y \) is uncorrelated with \( (X,C) \), so \( \beta_{X\mid C}=\theta+h\kappa \), where \( \kappa \) is the
coefficient of \( X \) in \( L(U\mid X,C) \). Put \( E=C-\delta X=gU+U_C \), uncorrelated with \( X \). As in the proof of
@thm-cau-bad-controls(b), \( L(U\mid X,C)=L(U\mid E)=\lambda(C-\delta X) \) with \( \lambda=g\sigma_U^2/(g^2\sigma_U^2+\sigma_C^2) \),
so \( \kappa=-\lambda\delta \). The mediator is a collider on the path \( X\to C\leftarrow U\to Y \), which adjustment opens.
:::

::: {#exr-cau-bias-amplification}
[B3]

Let \( U \) be unrecorded, \( X=\pi I+gU+U_X \) and \( Y=\tau X+hU+U_Y \), with \( I \), \( U \), \( U_X \), \( U_Y \) independent, \( \Var I=\Var U=1 \)
and \( \Var U_X=\sigma^2 \). Show that the slope of \( L(Y\mid X) \) is \( \tau+gh/(\pi^2+g^2+\sigma^2) \), and that the coefficient of \( X \) in
\( L(Y\mid X,I) \) is \( \tau+gh/(g^2+\sigma^2) \). Adjusting for the instrument \( I \) enlarges the bias.
:::

::: {.solution}
\( \Cov(X,Y)=\tau\Var X+gh \) and \( \Var X=\pi^2+g^2+\sigma^2 \), which gives the first slope. For the second, by
@prp-cor-partial-coefficient the coefficient is \( \sigma_{YX\cdot I}/\sigma_{XX\cdot I} \). Removing the linear effect of \( I \) leaves
\( X-\pi I=gU+U_X \), with variance \( g^2+\sigma^2 \), and \( Y-L(Y\mid I)=\tau(gU+U_X)+hU+U_Y \), so
\( \sigma_{YX\cdot I}=\tau(g^2+\sigma^2)+gh \). Divide. The confounded covariance \( gh \) is the same in both, but adjustment divides it by
a smaller variance.
:::

### C. Going deeper

::: {#exr-cau-butterfly}
[C1]

*(Butterfly bias.)* Let \( H_1,H_2,U_C,U_X,U_Y \) be independent with mean zero and variance \( 1 \), and let
\[
C=H_1+H_2+U_C,\qquad X=H_1+C+U_X,\qquad Y=\tau X+H_2+C+U_Y ,
\]
so that \( C \) is at once a collider on \( H_1\to C\leftarrow H_2 \) and a common cause of \( X \) and
\( Y \).

::: {.enumerate options="label=(\alph*)"}
1. List the back-door paths from \( X \) to \( Y \) and say which are open when \( C \) is not conditioned on,
   and which when it is. Conclude that neither \( \emptyset \) nor \( \{C\} \) satisfies the back-door
   criterion (@def-cau-backdoor-criterion).

2. Show that the slope of \( L(Y\mid X) \) is \( \tau+\tfrac57 \) and that the coefficient of \( X \) in
   \( L(Y\mid X,C) \) is \( \tau-\tfrac15 \). Adjusting for \( C \) reduces the bias and reverses its sign, so
   the two estimates bracket \( \tau \).

3. Now change one arrow, replacing \( C=H_1+H_2+U_C \) by \( C=H_1-H_2+U_C \) and keeping everything else.
   Show that the two biases become \( \tfrac37 \) and \( \tfrac15 \): both positive, so the two estimates no
   longer bracket the truth.
:::

There is therefore no "adjust and compare" rule that turns two biased estimates into bounds.
:::

::: {.solution}
(a) Write the paths that leave \( X \) by an arrow into \( X \). They are \( X\leftarrow C\to Y \),
\( X\leftarrow H_1\to C\to Y \), \( X\leftarrow C\leftarrow H_2\to Y \) and
\( X\leftarrow H_1\to C\leftarrow H_2\to Y \). On the first three, \( C \) is a non-collider, so they are open
when \( C \) is free and blocked when it is conditioned on. On the fourth, \( C \) is a collider, so that path
is blocked when \( C \) is free and opened when it is conditioned on. Condition (i) is not at issue, since
\( C \) is no descendant of \( X \); it is (ii) that fails for both candidate sets, because every choice
leaves one of the four paths open.

(b) Substituting, \( X=2H_1+H_2+U_C+U_X \) and \( Y=\tau X+R \) with \( R=H_1+2H_2+U_C+U_Y \). All five
disturbances are uncorrelated with variance one, so
\( \Var(X)=4+1+1+1=7 \), \( \Cov(X,R)=2+2+1=5 \), \( \Var(C)=3 \), \( \Cov(X,C)=2+1+1=4 \) and
\( \Cov(R,C)=1+2+1=4 \). The slope of \( L(Y\mid X) \) is \( \tau+\Cov(X,R)/\Var(X)=\tau+\tfrac57 \)
(@thm-proj-blp). For the adjusted coefficient, \( U_Y \) is uncorrelated with \( (X,C) \), so the coefficient
of \( X \) in \( L(Y\mid X,C) \) is \( \tau+k_1 \), where \( (k_1,k_2) \) solves the population normal
equations
\[
\begin{pmatrix}7&4\\4&3\end{pmatrix}\begin{pmatrix}k_1\\k_2\end{pmatrix}
=\begin{pmatrix}5\\4\end{pmatrix}.
\]
The determinant is \( 5 \) and \( k_1=(5\cdot3-4\cdot4)/5=-\tfrac15 \).

(c) Now \( X=2H_1-H_2+U_C+U_X \) with \( \Var(X)=7 \) again, and \( R=H_2+C+U_Y=H_1+U_C+U_Y \), so
\( \Cov(X,R)=2+1=3 \) and the unadjusted bias is \( \tfrac37 \). Also \( \Var(C)=3 \),
\( \Cov(X,C)=2+1+1=4 \) and \( \Cov(R,C)=1+1=2 \), so the same matrix now faces the right-hand side
\( (3,2)\T \) and \( k_1=(3\cdot3-4\cdot2)/5=\tfrac15 \).

The sign of the M-bias that conditioning opens is the sign of a product of four path coefficients
(@eq-cau-bias-m), while the sign of the confounding that conditioning closes is set by the arrows out of
\( C \); changing one arrow changes the first and not the second. Bracketing in (b) is an accident of the
numbers, and nothing in the graph guarantees it.
:::

::: {#exr-cau-cancelling-control}
[C2]

*(A bad control that does no harm.)* Let \( H \) be unrecorded, \( U_X,U_C,U_Y \) independent of it and of
each other, \( \Var(H)=\Var(U_X)=1 \), \( \Var(U_C)=\sigma_C^2>0 \), and
\[
X=H+U_X,\qquad Y=\tau X+H+U_Y,\qquad C=pH+qX+U_C .
\]
The control \( C \) is measured after the treatment: it is a descendant of \( X \) and a collider on
\( X\to C\leftarrow H\to Y \), so the advice of this section forbids it twice over.

::: {.enumerate options="label=(\alph*)"}
1. Show that the slope of \( L(Y\mid X) \) is \( \tau+\tfrac12 \).

2. Show that the coefficient of \( X \) in \( L(Y\mid X,C) \) is
   \[
   \tau+\frac{\sigma_C^2-pq}{p^2+2\sigma_C^2}.
   \]

3. Deduce that the adjusted coefficient is exactly \( \tau \) when \( \sigma_C^2=pq \), for instance at
   \( p=q=\sigma_C^2=1 \), and check the two extreme cases \( q=0 \) and \( p=0 \) against
   @exr-cau-descendant-only and against the idea of a proxy for \( H \).
:::

Why does this not overturn the advice?
:::

::: {.solution}
(a) \( Y=\tau X+R \) with \( R=H+U_Y \), and \( \Var(X)=2 \), \( \Cov(X,R)=1 \), so the slope is
\( \tau+\tfrac12 \).

(b) Since \( U_Y \) is uncorrelated with \( (X,C) \), the coefficient of \( X \) in \( L(Y\mid X,C) \) is
\( \tau+k_1 \), with \( (k_1,k_2) \) the coefficients of \( L(H\mid X,C) \). Write
\( C=(p+q)H+qU_X+U_C \). Then \( \Var(X)=2 \), \( \Cov(X,C)=p+2q \),
\( \Var(C)=(p+q)^2+q^2+\sigma_C^2 \), \( \Cov(H,X)=1 \) and \( \Cov(H,C)=p+q \). The determinant of the
covariance matrix of \( (X,C) \) is
\[
D=2\bigl\{(p+q)^2+q^2+\sigma_C^2\bigr\}-(p+2q)^2=p^2+2\sigma_C^2>0 ,
\]
and solving the two equations for \( k_1 \) gives
\( k_1=\bigl\{(p+q)^2+q^2+\sigma_C^2-(p+2q)(p+q)\bigr\}/D=(\sigma_C^2-pq)/D \), since the brace equals
\( (p+q)\{(p+q)-(p+2q)\}+q^2+\sigma_C^2=\sigma_C^2-pq \).

(c) At \( p=q=\sigma_C^2=1 \) the numerator vanishes and the adjusted coefficient is exactly \( \tau \), while
the unadjusted one is \( \tau+\tfrac12 \). At \( q=0 \), \( C=pH+U_C \) is a pre-treatment proxy for the
confounder and the bias is \( \sigma_C^2/(p^2+2\sigma_C^2)<\tfrac12 \), tending to \( 0 \) as
\( \sigma_C^2\to0 \): a perfect proxy is as good as \( H \) itself. At \( p=0 \), \( C=qX+U_C \) is a pure
descendant of \( X \) and the bias is \( \sigma_C^2/(2\sigma_C^2)=\tfrac12 \), unchanged, which is
@exr-cau-descendant-only again.

The advice survives because the cancellation is not a property of the graph. Conditioning on \( C \) opens
the path \( X\to C\leftarrow H\to Y \), whose contribution is \( -pq/D \), and shrinks the confounding along
\( X\leftarrow H\to Y \) to \( \sigma_C^2/D \); the two are equal only on the surface \( \sigma_C^2=pq \) in
the parameter space, and \( p \), \( q \) and \( \sigma_C^2 \) all involve the unrecorded \( H \), so no data
tell the analyst where on that surface the truth lies. A rule that forbids the control is right in the sense
that matters: it is right for every graph, while the exception holds for a null set of the parameters and
cannot be recognized.
:::
