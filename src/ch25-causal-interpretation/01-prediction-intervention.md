# Prediction versus intervention

A regression of \( Y \) on \( X \) answers a question about *seeing*: among units in which \( X=x \) was
observed, what is \( Y \) typically? The answer is \( \E(Y\mid X=x) \), or its best linear approximation
\( L(Y\mid X) \) (@thm-proj-blp), which least squares estimates whatever produced the data (@prp-proj-consistency). A decision asks about *doing*: if \( X \) were set to \( x \) by an outside agent, what
would \( Y \) typically be? This section makes the second question precise and shows, in the simplest model,
that the answers differ, that more data do not close the gap, and that the better predictor can be the worse
guide to action.

## Structural causal models

A joint distribution says how variables co-vary, not which would move if another were changed. For that we need
a model of the mechanism, a list of the processes that generate each variable. The following definition goes back
to Wright (1921) and Haavelmo (1943).

::: {#def-cau-scm}
[Structural causal model]

Let \( V_1,\dots,V_m \) be random variables. A **structural causal model** (SCM) for them consists of

::: {.enumerate options="label=(\alph*)"}
1. for each \( j \), a set \( \mathrm{pa}(j)\subseteq\{1,\dots,j-1\} \) of indices, the **parents** of \( V_j \);

2. for each \( j \), a function \( f_j \), and mutually independent random variables \( U_1,\dots,U_m \), the
   **disturbances**, such that
   \[
   V_j=f_j\bigl(V_{\mathrm{pa}(j)},U_j\bigr),\qquad j=1,\dots,m,
   \]{#eq-cau-scm}

   where \( V_{\mathrm{pa}(j)} \) is the vector of parents of \( V_j \).
:::

Each line of @eq-cau-scm is an **assignment**: it says how nature sets \( V_j \), not merely that an equation
holds. Solving the assignments in the order \( j=1,\dots,m \) writes every \( V_j \) as a function of
\( U_1,\dots,U_j \). For a variable \( X=V_k \) and a value \( x \), the **intervention**
\( \operatorname{do}(X=x) \) replaces the \( k \)th assignment by \( V_k=x \) and keeps every other assignment and
every disturbance unchanged. The distribution of the variables in the modified model is the
**interventional distribution**, and expectations under it are written \( \E\bigl(Y\mid\operatorname{do}(X=x)\bigr) \).
The model is **linear** if each \( f_j(\bv,u)=c_j+\sum_{l\in\mathrm{pa}(j)}b_{jl}v_l+u \) with constants
\( c_j,b_{jl} \), and disturbances with mean zero and finite variance.
:::

[Section 25.3](03-confounding-dags.html) draws the parent sets as a graph. The substantive assumption is that an
intervention changes one assignment and leaves the other mechanisms intact. Conditioning and intervening are then
different operations: \( \E(Y\mid X=x) \) looks at the part of the *original* distribution where \( X=x \) happened,
while \( \E\bigl(Y\mid\operatorname{do}(X=x)\bigr) \) describes a *new* distribution in which \( X=x \) was imposed.

## A common cause

The simplest model in which the two differ has a background variable \( Z \) that influences both \( X \) and \( Y \):
\[
Z=c_Z+U_Z,\qquad X=c_X+\alpha Z+U_X,\qquad Y=c_Y+\theta X+\gamma Z+U_Y,
\]{#eq-cau-confounder-model}

with independent disturbances of variances \( \sigma_Z^2>0 \), \( \sigma_X^2>0 \) and \( \sigma_Y^2 \)
([Figure 25.1.1](#fig-cau-confounder-dag)). The coefficient \( \theta \) is what the assignment for \( Y \) does with
its argument \( X \).

\begin{center}
\begin{tikzpicture}[>=Stealth, line cap=round]
\definecolor{ink}{HTML}{1B1F27}\definecolor{accent}{HTML}{2F5F96}\definecolor{draft}{HTML}{A45C25}
\tikzset{var/.style={draw=ink, circle, minimum size=7.5mm, inner sep=0pt, font=\small}}
\node[var, dashed] (Z) at (1.6,1.25) {$Z$};
\node[var] (X) at (0,0) {$X$};
\node[var] (Y) at (3.2,0) {$Y$};
\draw[->, thick, draft] (Z) -- node[above left, font=\small] {$\alpha$} (X);
\draw[->, thick, draft] (Z) -- node[above right, font=\small] {$\gamma$} (Y);
\draw[->, thick, accent] (X) -- node[below, font=\small] {$\theta$} (Y);
\end{tikzpicture}
\end{center}

[**Figure 25.1.1.** The model @eq-cau-confounder-model. An arrow points from each parent to its child and is
labelled with its coefficient. The dashed circle marks a variable that is not recorded in @exm-cau-visits.]{#fig-cau-confounder-dag}

::: {#prp-cau-confounded}
[Seeing and doing with a common cause]

In the model @eq-cau-confounder-model:

::: {.enumerate options="label=(\alph*)"}
1. \( \E\bigl(Y\mid\operatorname{do}(X=x)\bigr)=c_Y+\gamma c_Z+\theta x \), so the interventional slope is \( \theta \);

2. the best linear predictor \( L(Y\mid X) \) has slope
   \[
   \beta^*=\theta+\frac{\gamma\alpha\sigma_Z^2}{\alpha^2\sigma_Z^2+\sigma_X^2};
   \]{#eq-cau-confounded-slope}

3. the coefficient of \( X \) in \( L(Y\mid X,Z) \) is \( \theta \);

4. if the disturbances are normal, \( \E(Y\mid X=x)=L(Y\mid X) \) evaluated at \( x \), so the slopes of seeing and
   doing differ whenever \( \gamma\alpha\ne0 \).
:::

:::

::: {.proof}
(a) Under \( \operatorname{do}(X=x) \) the assignment for \( Z \) is unchanged and does not involve \( X \), so \( Z=c_Z+U_Z \)
still, and \( Y=c_Y+\theta x+\gamma Z+U_Y \). Take expectations.
(b) By @thm-proj-blp the slope is \( \Cov(X,Y)/\Var(X) \). Here \( \Var(X)=\alpha^2\sigma_Z^2+\sigma_X^2 \) and
\( \Cov(X,Y)=\theta\Var(X)+\gamma\Cov(X,Z)=\theta\Var(X)+\gamma\alpha\sigma_Z^2 \).
(c) The error \( Y-c_Y-\theta X-\gamma Z=U_Y \) has mean zero and is uncorrelated with \( X \) and \( Z \), which are
functions of \( U_Z \) and \( U_X \). So \( c_Y+\theta X+\gamma Z \) satisfies the population normal equations of
@thm-proj-blp and is the projection.
(d) \( (X,Y) \) is a linear function of the normal vector \( (U_Z,U_X,U_Y) \), hence bivariate normal, and its
conditional mean is the best linear predictor (@thm-mvn-conditional).
:::

The extra term in @eq-cau-confounded-slope is association that \( X \) borrows from \( Z \). Part (c) is the familiar
remedy, including the common cause as a regressor, and it needs \( Z \) to be recorded.
[Section 25.4](04-backdoor.html) says in general which variables must be included.

::: {#exm-cau-visits}
[Clinic visits and health]

In a synthetic population, \( Z \) is the severity of a chronic illness, which nobody records; \( X \) is the
number of clinic visits in a year; and \( Y \) is a health score measured at the end of the year. Sicker
people visit more (\( \alpha=1.5 \)) and end the year less healthy (\( \gamma=-6 \)), and each visit helps
(\( \theta=2 \)). With \( \sigma_Z=\sigma_X=1 \) and \( \sigma_Y=2 \), @eq-cau-confounded-slope gives
\[
\beta^*=2+\frac{(-6)(1.5)}{1.5^2+1}=2-\frac{9}{3.25}=-0.769 .
\]
In a sample of \( n=2000 \) people the least squares slope is \( -0.790 \): each extra visit
goes with a *lower* health score. A second data set in which visits are assigned at random, which is an
intervention on \( X \), has least squares slope \( 1.999 \)
([Figure 25.1.2](#fig-cau-visits)). A health authority that capped visits on the strength of the first regression
would harm its patients.
:::

::: {when-format="html"}
![**Figure 25.1.2.** The model of @exm-cau-visits. (a) Observational data: the fitted line slopes down. (b) Data in
which the number of visits is set at random by the analyst: the same mechanism for \( Y \), but now \( X \) carries no
information about illness. The solid line is the observational fit in both panels; the dashed line has the causal
slope \( \theta=2 \).](prediction_intervention.svg){#fig-cau-visits width=100%}
:::

::: {when-format="pdf"}
![The model of @exm-cau-visits. (a) Observational data: the fitted line slopes down. (b) Data in
which the number of visits is set at random by the analyst: the same mechanism for \( Y \), but now \( X \) carries no
information about illness. The solid line is the observational fit in both panels; the dashed line has the causal
slope \( \theta=2 \).](prediction_intervention.pdf){width=100%}
:::

```{.python .run #cell-prediction-simulate}
import numpy as np

rng = np.random.default_rng(2501)
n = 2000
alpha, theta, gamma = 1.5, 2.0, -6.0          # Z -> X, X -> Y, Z -> Y

def draw(n, do_x=None):
    """Draw from the model; with do_x, X is set by the analyst instead of by Z."""
    Z = rng.normal(size=n)                     # illness severity, never recorded
    X = 4 + alpha * Z + rng.normal(size=n) if do_x is None else do_x
    Y = 50 + theta * X + gamma * Z + rng.normal(scale=2.0, size=n)
    return X, Y

x_obs, y_obs = draw(n)                                  # observational data
x_do = rng.uniform(x_obs.min(), x_obs.max(), size=n)    # X assigned at random
x_int, y_int = draw(n, do_x=x_do)                       # interventional data

slope_obs, icpt_obs = np.polyfit(x_obs, y_obs, 1)
slope_int, icpt_int = np.polyfit(x_int, y_int, 1)
print(f"slope from observational data:  {slope_obs:.3f}")
print(f"slope from interventional data: {slope_int:.3f}   (theta = {theta})")
```

## A good predictor can be a bad causal model

The observational line in @exm-cau-visits is the best linear predictor of health from visits in the population
that generated the data, and the causal line, with slope \( \theta \) through the same means, predicts *worse*
there. On the observational data the fitted line has mean
squared error \( 15.0 \) and the causal line \( 40.7 \), which is worse than predicting every
person by the sample mean (the variance of \( Y \) is \( 17.1 \)). On the interventional data the
positions reverse: \( 123.3 \) for the fitted line and \( 39.2 \) for the causal line.

```{.python .run #cell-prediction-mse}
def mse(x, y, slope, icpt):
    return np.mean((y - icpt - slope * x) ** 2)

causal_icpt = np.mean(y_obs) - theta * np.mean(x_obs)   # causal slope, fitted level
for label, x, y in [("observational", x_obs, y_obs), ("interventional", x_int, y_int)]:
    print(f"{label:15s} MSE: regression line {mse(x, y, slope_obs, icpt_obs):6.2f},"
          f" causal line {mse(x, y, theta, causal_icpt):6.2f}")
```

Observationally, any line through the means other than the projection loses exactly
\( (\text{slope}-\beta^*)^2\Var(X) \) (@exr-cau-excess-mse). Under intervention a large \( X \) no longer signals a large
\( Z \), so the part of the slope borrowed from \( Z \) becomes pure error: the causal line describes the mechanism both regimes
share, the regression line a distribution that the intervention destroys.

Confounding by a common cause is one of three standard ways in which a predictor misleads about action. In
**reverse causation** \( Y \) causes \( X \), and \( X \) can predict \( Y \) perfectly while setting it does nothing (@exr-cau-reverse): a thermometer predicts the temperature, but moving its needle does not warm the room. In
**selection** the data contain only units selected by a common effect of \( X \) and \( Y \), which creates association
without effect ([Section 25.5](05-bad-controls.html)).

More data do not help: the least squares slope converges to \( \beta^* \), not \( \theta \) (@prp-proj-consistency), and
no procedure using only the distribution of \( (X,Y) \) can do better, because models with different \( \theta \) can produce the same distribution (@exr-cau-observational-equivalence). Causal conclusions need either an intervention, as in panel (b) of
[Figure 25.1.2](#fig-cau-visits), or assumptions about the mechanism that link the interventional distribution to the
observational one.

## When seeing is doing

Conditioning and intervening agree for every mechanism when nothing in the model influences \( X \).

::: {#prp-cau-exogenous}
[An exogenous cause]

In a structural causal model, suppose \( X=V_k \) has no parents, so that \( X=f_k(U_k) \). Then a version of the
conditional distribution of the other variables given \( X=x \) is their distribution under \( \operatorname{do}(X=x) \); in
particular, when \( \E|Y|<\infty \), \( \E(Y\mid X=x)=\E\bigl(Y\mid\operatorname{do}(X=x)\bigr) \) for almost every \( x \) (with respect to the distribution
of \( X \)), and for every \( x \) with \( \Pr(X=x)>0 \).
:::

::: {.proof}
Solve the assignments in order. Every variable is a function \( h(X,\mathbf{U}_{-k}) \) of \( X \) and the disturbances
\( \mathbf{U}_{-k} \) other than \( U_k \), because \( U_k \) enters the other assignments only through \( X \). The same
function describes the model under \( \operatorname{do}(X=x) \), with \( x \) in place of \( X \), since only the
\( k \)th assignment changed. Now \( X=f_k(U_k) \) is independent of \( \mathbf{U}_{-k} \), and for independent \( X \)
and \( \mathbf{U}_{-k} \) a version of the conditional distribution of \( h(X,\mathbf{U}_{-k}) \) given \( X=x \) is the distribution of
\( h(x,\mathbf{U}_{-k}) \) (Billingsley 1995), which is the interventional distribution.
:::

A randomized experiment makes the treatment exogenous by construction: the random device is the only parent of \( X \)
and influences nothing else, so what is seen is what would be done. [Section 25.2](02-potential-outcomes.html) turns
this into exact finite-sample statements that need no structural model.

## Exercises

### A. Check your understanding

::: {#exr-cau-zero-slope}
[A1]

In the model @eq-cau-confounder-model with \( \theta=2 \), \( \gamma=-6 \) and \( \sigma_Z=\sigma_X=1 \), for which
values of \( \alpha \) is the observational slope \( \beta^* \) exactly zero? What would an analyst who saw
\( \beta^*=0 \) be tempted to conclude?
:::

### B. Practice

::: {#exr-cau-excess-mse}
[B1]

Let \( (x_i,y_i) \), \( i=1,\dots,n \), be any data with least squares slope \( \hat{\beta} \). For a number \( b \),
let \( \ell_b(x)=\bar{y}+b(x-\bar{x}) \) be the line with slope \( b \) through the means. Show that
\[
\frac1n\sum_i\bigl(y_i-\ell_b(x_i)\bigr)^2=\frac1n\sum_i\bigl(y_i-\ell_{\hat{\beta}}(x_i)\bigr)^2
+(b-\hat{\beta})^2\,\frac1n\sum_i(x_i-\bar{x})^2 ,
\]
and deduce the ordering of the two observational mean squared errors in @exm-cau-visits.
:::

::: {.solution}
Write \( y_i-\ell_b(x_i)=\bigl(y_i-\ell_{\hat{\beta}}(x_i)\bigr)+(\hat{\beta}-b)(x_i-\bar{x}) \). The cross term is
\( (\hat{\beta}-b)\sum_i\hat{e}_i(x_i-\bar{x}) \), where \( \hat{e}_i \) are the least squares residuals. The residuals are
orthogonal to \( \bone \) and to \( \x \) (@thm-proj-ls-projection), so the cross term is zero, and squaring gives the identity.
With \( b=\theta=2 \) the causal line through the means loses \( (2-\hat{\beta})^2 \) times the variance of the visits,
which is why its mean squared error exceeds that of the fitted line. The script checks the identity to rounding error.
:::

::: {#exr-cau-reverse}
[B2]

Let \( Y=c+U_Y \) and \( X=\delta Y+U_X \), with independent disturbances of variances \( \sigma_Y^2 \) and
\( \sigma_X^2 \). Find the slope of \( L(Y\mid X) \), the squared correlation of \( X \) and \( Y \), and
\( \E\bigl(Y\mid\operatorname{do}(X=x)\bigr) \). Show that the squared correlation can be arbitrarily close to one while
setting \( X \) has no effect on \( Y \).
:::

::: {#exr-cau-observational-equivalence}
[B3]

Take \( c_Z=c_X=c_Y=0 \) in @eq-cau-confounder-model with normal disturbances, and suppose only \( (X,Y) \) is observed.
Show that the parameter values
\[
(\alpha,\theta,\gamma,\sigma_Z^2,\sigma_X^2,\sigma_Y^2)=(1,1,0,1,1,1)
\quad\text{and}\quad
\bigl(1,0,\tfrac43,\tfrac32,\tfrac12,\tfrac13\bigr)
\]
give the same distribution of \( (X,Y) \), although \( X \) has an effect on \( Y \) in the first model and none in the
second.
:::

::: {.solution}
Both distributions are bivariate normal with mean zero, so it suffices to compare covariance matrices, using
\( \Var(X)=\alpha^2\sigma_Z^2+\sigma_X^2 \), \( \Cov(X,Y)=\theta\Var(X)+\gamma\alpha\sigma_Z^2 \) and
\( \Var(Y)=\theta^2\Var(X)+\gamma^2\sigma_Z^2+2\theta\gamma\alpha\sigma_Z^2+\sigma_Y^2 \). In the first model these
are \( 2 \), \( 2 \) and \( 2+1=3 \). In the second they are \( \tfrac32+\tfrac12=2 \), \( \tfrac43\cdot\tfrac32=2 \) and
\( \tfrac{16}9\cdot\tfrac32+\tfrac13=3 \). The covariance matrices agree, so no function of the \( (X,Y) \) data can
distinguish \( \theta=1 \) from \( \theta=0 \).
:::

### C. Going deeper

::: {#exr-cau-any-effect}
[C1]

Take \( c_Z=c_X=c_Y=0 \) and normal disturbances in @eq-cau-confounder-model, and let \( \bSigma \) be any
\( 2\times2 \) positive definite matrix, with entries \( v=\bSigma_{11} \), \( c=\bSigma_{12} \) and
\( w=\bSigma_{22} \). Show that for **every** real \( \theta \) there are values of
\( (\alpha,\gamma,\sigma_Z^2,\sigma_X^2,\sigma_Y^2) \), with \( \sigma_Z^2>0 \), \( \sigma_X^2>0 \) and
\( \sigma_Y^2\ge0 \), for which \( (X,Y)\sim\Normal_2(\bzero,\bSigma) \). Hence the observational distribution
constrains the causal effect not at all: @exr-cau-observational-equivalence exhibits two values of \( \theta \)
compatible with one distribution, and in fact no value is ever excluded.
:::

::: {.solution}
Put \( \sigma_Z^2=1 \) and \( \alpha=\sqrt t \) for a number \( t\in(0,v) \) to be chosen, and
\( \sigma_X^2=v-t>0 \), so that \( \Var(X)=\alpha^2\sigma_Z^2+\sigma_X^2=v \) whatever \( t \) is. Take
\( \gamma=(c-\theta v)/\sqrt t \), so that, by the covariance formulas in the solution to
@exr-cau-observational-equivalence,
\( \Cov(X,Y)=\theta v+\gamma\alpha\sigma_Z^2=\theta v+(c-\theta v)=c \). The same formulas leave
\[
\sigma_Y^2=w-\theta^2v-\gamma^2-2\theta(c-\theta v)
=w+\theta^2v-2\theta c-\frac{(c-\theta v)^2}{t}=:\varphi(t).
\]
As \( t\uparrow v \), and using \( (c-\theta v)^2/v=c^2/v-2\theta c+\theta^2v \),
\[
\varphi(t)\ \longrightarrow\ w+\theta^2v-2\theta c-\frac{c^2}{v}+2\theta c-\theta^2v=w-\frac{c^2}{v}>0 ,
\]
the last inequality because \( \bSigma \) is positive definite. So \( \varphi(t)>0 \) for some \( t<v \); fix
such a \( t \) and set \( \sigma_Y^2=\varphi(t) \). All three variances are admissible, and \( (X,Y) \) is a
linear function of a normal vector, hence normal (@thm-mvn-linear) with mean zero and the prescribed
covariance matrix. Every \( \theta \) is therefore compatible with every observational distribution: not only
can the sign of the effect not be recovered, neither can its size.
:::

::: {#exr-cau-direction-identified}
[C2]

The impossibility in @exr-cau-any-effect is a fact about normal disturbances, not about observational data.
Let \( U_X \) and \( U_Y \) be independent, each taking the values \( \pm1 \) with probability \( \tfrac12 \), and
let the mechanism be \( X=U_X \), \( Y=X+U_Y \). Show that no linear structural causal
model (@def-cau-scm) in which \( Y \) has no parents and \( X=c+\delta Y+U \), with \( U \) independent of \( Y \),
reproduces the joint distribution of \( (X,Y) \), whatever the law of \( U \) and whatever \( c \) and
\( \delta \) are. So here the direction of the arrow *is* visible in the data.
:::

::: {.solution}
Under the stated mechanism \( (X,Y) \) takes the four values \( (1,2) \), \( (1,0) \), \( (-1,0) \), \( (-1,-2) \),
each with probability \( \tfrac14 \); in particular \( Y\in\{-2,0,2\} \) with probabilities
\( \tfrac14,\tfrac12,\tfrac14 \). Suppose \( X=c+\delta Y+U \) with \( U \) independent of \( Y \). Then the
conditional law of \( X-c-\delta Y \) given \( Y=y \) is the law of \( U \) for every \( y \) with
\( \Pr(Y=y)>0 \). But given \( Y=2 \) the only possible value of \( X \) is \( 1 \), so that conditional law is
degenerate at \( 1-c-2\delta \); while given \( Y=0 \), \( X \) equals \( 1 \) or \( -1 \) with probability
\( \tfrac12 \) each, so the conditional law puts mass \( \tfrac12 \) on each of \( 1-c \) and \( -1-c \) and is
not degenerate. A degenerate law and a two-point law cannot both be the law of \( U \), a contradiction.

The contrast with @exr-cau-any-effect is the whole point. Two independent non-normal disturbances leave a
trace in the joint law that tells a cause from an effect; two independent normal ones do not, because a
normal vector is determined by its first two moments and those are symmetric in the two directions. Methods
that exploit non-normality to orient arrows rest on this observation.
:::

::: {#exr-cau-invariance}
[C3]

In @eq-cau-confounder-model, suppose an outside agent sets \( X \) to a random value \( \tilde X \) with mean
\( m \) and variance \( \nu>0 \), independent of \( (U_Z,U_X,U_Y) \); this is
\( \operatorname{do}(X=\tilde X) \), and the observational regime is *not* of this form.

::: {.enumerate options="label=(\alph*)"}
1. Show that in the new regime the line \( \ell(x)=a+bx \) has mean squared prediction error
   \[
   \E\bigl\{Y-\ell(\tilde X)\bigr\}^2
   =\bigl\{k+(\theta-b)m\bigr\}^2+(\theta-b)^2\nu+\gamma^2\sigma_Z^2+\sigma_Y^2 ,
   \qquad k=c_Y+\gamma c_Z-a .
   \]

2. Deduce that the minimum is attained exactly at \( b=\theta \) and \( a=c_Y+\gamma c_Z \), whatever the
   agent's device is, and that the minimum value \( \gamma^2\sigma_Z^2+\sigma_Y^2 \) does not depend on the
   device either.

3. Deduce that the observational line of @eq-cau-confounded-slope has excess error at least
   \( (\beta^*-\theta)^2\nu \), which grows without bound as the agent spreads \( X \) more widely, and that a
   line that is optimal under two devices with different means is the causal line.
:::

The causal line is thus the only line that predicts well under every intervention; the regression line is
optimal for one regime, the one the intervention destroys.
:::

::: {.solution}
(a) Under \( \operatorname{do}(X=\tilde X) \) the assignments for \( Z \) and \( Y \) are unchanged, so
\( Z=c_Z+U_Z \) and \( Y=c_Y+\theta\tilde X+\gamma Z+U_Y \). Hence
\[
Y-\ell(\tilde X)=k+(\theta-b)\tilde X+\gamma U_Z+U_Y ,
\]
with \( k=c_Y+\gamma c_Z-a \). The last two terms have mean zero and are independent of \( \tilde X \), so the
mean squared error is \( \E\{k+(\theta-b)\tilde X\}^2+\gamma^2\sigma_Z^2+\sigma_Y^2 \), and
\( \E\{k+(\theta-b)\tilde X\}^2=\{k+(\theta-b)m\}^2+(\theta-b)^2\nu \).

(b) Both squared terms are nonnegative and \( \nu>0 \), so the second forces \( b=\theta \) and the first then
forces \( k=0 \), that is \( a=c_Y+\gamma c_Z \). Neither the minimizer nor the minimum involves \( m \) or
\( \nu \).

(c) With \( b=\beta^* \) the second term alone is \( (\beta^*-\theta)^2\nu \), and by
@eq-cau-confounded-slope this is nonzero whenever \( \gamma\alpha\ne0 \); it grows linearly in \( \nu \). If a
line is optimal for devices with means \( m_1\ne m_2 \) and positive variances, (b) applied to each gives
\( b=\theta \) and \( a=c_Y+\gamma c_Z \) — the same line. The observational regime escapes (b) only because
there \( X \) is *not* independent of \( U_Z \): the association it borrows from \( Z \) is real, and worth
using, exactly as long as nobody intervenes.
:::
