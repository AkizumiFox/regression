# Adjustment sets and the back-door criterion

If confounding is association along back-door paths, the remedy is to block them by conditioning on variables on
them, that is, by including those variables as regressors. This section says which sets do the job, proves that they do,
and compares their precision.

## The criterion

::: {#def-cau-backdoor-criterion}
[Back-door criterion]

Let \( X \) and \( Y \) be nodes of a DAG. A set \( \mathcal Z \) of nodes, containing neither \( X \) nor \( Y \), satisfies the
**back-door criterion** relative to \( (X,Y) \) if

::: {.enumerate options="label=(\roman*)"}
1. no node in \( \mathcal Z \) is a descendant of \( X \), and

2. \( \mathcal Z \) blocks every back-door path from \( X \) to \( Y \) (@def-cau-confounding).
:::

Such a set is called a (back-door) **adjustment set**.
:::

Condition (ii) closes the non-causal routes. Condition (i) protects the causal ones: a descendant of \( X \) may be a
mediator, whose conditioning removes part of the effect, or a collider, whose conditioning opens a path
([Section 25.5](05-bad-controls.html)). The empty set qualifies exactly when the effect is not confounded, and the parents
of \( X \) always qualify, provided \( Y \) is not one of them, since every back-door path begins \( X\leftarrow P \) with \( P \) a
parent and a non-collider. The
graph itself is an assumption that the data cannot check (@exr-cau-observational-equivalence); the criterion turns it into
a precise statement about which regression to run.

Recall from [Section 25.2](02-potential-outcomes.html) that \( Y(x) \) is the value of \( Y \) in the model modified by
\( \operatorname{do}(X=x) \), and that \( Y=Y(x) \) on the event \( X=x \). Conditional independence of \( Y(x) \) and \( X \) given
\( \mathbf{Z} \) is the potential-outcome form of "no confounding given \( \mathbf{Z} \)", which Rosenbaum and Rubin (1983), adding
positivity, called *strong ignorability*.

::: {#thm-cau-backdoor}
[Back-door adjustment]

Let \( \mathbf{V} \) follow a structural causal model with graph \( G \), let \( X \) and \( Y \) be two of its variables, and let
\( \mathbf{Z} \) be the vector of the variables in a set \( \mathcal Z \) of other nodes.

::: {.enumerate options="label=(\alph*)"}
1. **Adjustment formula.** Suppose \( X \) takes finitely many values, \( \Pr(X=x\mid\mathbf{Z})>0 \) with probability one,
   \( \E|Y(x)|<\infty \), and \( Y(x) \) is conditionally independent of \( X \) given \( \mathbf{Z} \). Then
   \[
   \E\bigl(Y\mid\operatorname{do}(X=x)\bigr)=\E\bigl[\E(Y\mid X=x,\mathbf{Z})\bigr],
   \]{#eq-cau-adjustment}

   where the outer expectation is over the distribution of \( \mathbf{Z} \).

2. **The criterion gives ignorability.** If \( \mathcal Z \) satisfies the back-door criterion relative to \( (X,Y) \), then
   \( Y(x) \) is conditionally independent of \( X \) given \( \mathbf{Z} \), for every \( x \).

3. **Linear models.** In a linear structural causal model whose disturbances have finite positive variances, if
   \( \mathcal Z \) satisfies the back-door criterion, the coefficient of \( X \) in \( L(Y\mid X,\mathbf{Z}) \) is the total effect
   \( \tau \) of \( X \) on \( Y \). Consequently, for independent copies of \( \mathbf{V} \), the least squares coefficient of \( X \) in the
   regression of \( Y \) on \( 1 \), \( X \) and \( \mathbf{Z} \) converges to \( \tau \) with probability one.
:::

The proofs of (b) and (c) below are complete when \( \mathcal Z \) contains every parent of \( X \); for general adjustment
sets they are sketches, since the proof of (c) uses (b).
:::

::: {.proof}
(a) Using the tower rule, then the conditional independence together with positivity, then \( Y=Y(x) \) on \( \{X=x\} \),
\[
\begin{aligned}
\E\,Y(x)&=\E\bigl[\E(Y(x)\mid\mathbf{Z})\bigr]\\
&=\E\bigl[\E(Y(x)\mid X=x,\mathbf{Z})\bigr]=\E\bigl[\E(Y\mid X=x,\mathbf{Z})\bigr],
\end{aligned}
\]
and \( \E\,Y(x)=\E\bigl(Y\mid\operatorname{do}(X=x)\bigr) \) by definition.

(b) *When \( \mathcal Z \) contains the parents of \( X \).* Let \( X=V_k \), and write \( \mathbf{U}_{-k} \) for the disturbances
other than \( U_k \). A node that is not a descendant of \( X \) has no directed path from \( X \), so solving the assignments shows
that it is a function of \( \mathbf{U}_{-k} \) alone; by (i), so is \( \mathbf{Z} \). In the model modified by
\( \operatorname{do}(X=x) \), the disturbance \( U_k \) appears nowhere, so \( Y(x) \) is also a function of \( \mathbf{U}_{-k} \). Hence
\( U_k \) is independent of \( (Y(x),\mathbf{Z}) \), and given \( \mathbf{Z} \) it remains independent of \( Y(x) \), with its unconditional
distribution. Finally \( X=f_k(\mathbf{Z}_{\mathrm{pa}},U_k) \), where \( \mathbf{Z}_{\mathrm{pa}} \) is the part of \( \mathbf{Z} \) formed by the
parents of \( X \). Given \( \mathbf{Z} \), \( X \) is a function of \( U_k \) alone, and so is conditionally independent of \( Y(x) \).

*General adjustment sets (sketch).* Form a graph containing both the original variables and their counterparts in the
model modified by \( \operatorname{do}(X=x) \), joined through shared disturbance nodes (a *twin network*); \( Y(x) \) is a node of it,
and \( X \) has no arrows into the counterpart variables. A path from \( X \) to \( Y(x) \) that leaves \( X \) through a child of \( X \) must
turn back somewhere, at the latest at a disturbance node shared with a counterpart variable, so it meets a collider that is a
descendant of \( X \); by (i) neither that collider nor its descendants are in \( \mathcal Z \), and the path is blocked. A path that
leaves \( X \) through a parent starts like a back-door path, and one checks, using (ii) and again (i), that \( \mathcal Z \) blocks it.
So \( \mathcal Z \) d-separates \( X \) from \( Y(x) \), and @thm-cau-d-separation(a) gives the conditional independence. Pearl (1995, 2009) gives the criterion and its proof in
the language of interventions, and Shpitser, VanderWeele and Robins (2010) characterize exactly which sets satisfy
@eq-cau-adjustment for every model on \( G \).

(c) By @prp-cau-confounding-bias(a), \( Y=\tau X+R \), where \( R \) is the value of \( Y \) under \( \operatorname{do}(X=0) \), that is,
\( R=Y(0) \). The variables \( R \), \( X \) and \( \mathbf{Z} \) are linear functions of the disturbances with coefficients determined
by \( \B \), so their covariance matrix depends only on \( \B \) and \( \boldsymbol{\Omega} \). The covariance matrix of
\( (X,\mathbf{Z}) \) is positive definite, because \( \Cov(\mathbf{V})=\A\boldsymbol{\Omega}\A\T \) with \( \A \) invertible and
\( \boldsymbol{\Omega} \) positive definite. By @prp-cor-partial-coefficient the coefficient of \( X \) in \( L(Y\mid X,\mathbf{Z}) \) is
\[
\frac{\sigma_{YX\cdot Z}}{\sigma_{XX\cdot Z}}=\tau+\frac{\sigma_{RX\cdot Z}}{\sigma_{XX\cdot Z}} ,
\]
with \( \sigma_{XX\cdot Z}>0 \). Consider the linear model with the same \( \B \) and \( \boldsymbol{\Omega} \) but normal disturbances.
It has the same covariance matrix of \( (R,X,\mathbf{Z}) \), hence the same partial covariance \( \sigma_{RX\cdot Z} \). In it
\( (R,X,\mathbf{Z}) \) is jointly normal and, by (b), \( R=Y(0) \) is conditionally independent of \( X \) given \( \mathbf{Z} \), so
\( \sigma_{RX\cdot Z}=0 \) by @prp-mvn-partial-meaning(b). Therefore the coefficient is \( \tau \). (This step is complete
when \( \mathcal Z \) contains the parents of \( X \), and otherwise inherits the sketch of (b).) The last statement is @prp-proj-consistency.
:::

Part (a) averages the conditional mean at \( X=x \) over the population distribution of \( \mathbf{Z} \), not over its
distribution among units with \( X=x \) (direct standardization), and needs positivity: where \( X=x \) never occurs, the data
say nothing about \( \E(Y\mid X=x,\mathbf{Z}) \). Part (c) answers the question left open after the Frisch–Waugh–Lovell
theorem: the multiple-regression coefficient of \( X \) in a linear model is the effect of changing \( X \) whenever the other
regressors form an adjustment set.

::: {.warning}
Only the coefficient of \( X \) acquires a causal meaning. The coefficients of the adjustment variables are generally not
their effects, because \( X \) may mediate them; reporting them as effects is the "Table 2 fallacy" (Westreich and Greenland
2013), computed in @exr-cau-table-two.
:::

## An example with eight variables

::: {#exm-cau-eight}
[Choosing an adjustment set]

[Figure 25.4.1](#fig-cau-eight-dag) shows a linear model with unit-variance disturbances and
\[
\begin{aligned}
W_2&=0.8W_1+U,\qquad X=0.7W_1+0.6W_3+0.9I+U,\qquad M=0.8X+U,\\
Y&=0.5X+0.6M+0.7W_2+0.8W_3+P+U,
\end{aligned}
\]
with a separate disturbance \( U \) in each equation and \( W_1,W_3,I,P \) parentless. Summing over the two directed paths,
\( \tau=0.5+0.8\times0.6=0.98 \). The back-door paths \( X\leftarrow W_1\to W_2\to Y \) and \( X\leftarrow W_3\to Y \) are blocked by
\( W_1 \) or \( W_2 \), and by \( W_3 \), and \( M \) is a descendant of \( X \). So a set is valid iff it contains \( W_3 \) and one of
\( W_1,W_2 \), and not \( M \); \( I \) and \( P \) are optional.

The listing computes the coefficient of \( X \) in \( L(Y\mid X,\mathbf{Z}) \) exactly from \( \Cov(\mathbf{V})=\A\A\T \), with the
factor \( \sigma_{Y\cdot XZ}/\sigma_{X\cdot Z} \) that governs precision (see below), and the standard deviation over
\( 4000 \) simulated samples of size \( n=500 \):

| adjustment set | valid | coefficient of \( X \) | \( \sigma_{Y\cdot XZ}/\sigma_{X\cdot Z} \) | simulated sd |
|---|---|---|---|---|
| none | no | \( 1.308 \) | \( 1.15 \) | \( 0.052 \) |
| \( W_3 \) | no | \( 1.150 \) | \( 1.16 \) | \( 0.052 \) |
| \( W_1,W_3 \) | yes | \( 0.980 \) | \( 1.25 \) | \( 0.056 \) |
| \( W_2,W_3 \) | yes | \( 0.980 \) | \( 1.06 \) | \( 0.047 \) |
| \( W_2,W_3,P \) | yes | \( 0.980 \) | \( 0.80 \) | \( 0.036 \) |
| \( W_1,W_3,I \) | yes | \( 0.980 \) | \( 1.69 \) | \( 0.076 \) |
| \( W_1,W_3,M \) | no | \( 0.500 \) | \( 1.72 \) | \( 0.077 \) |
| \( W_1,W_2,W_3,I,P \) | yes | \( 0.980 \) | \( 1.17 \) | \( 0.052 \) |

Every valid set gives \( \tau \); the invalid sets leave a back-door path open or, with \( M \), keep only the direct effect
\( 0.5 \). Among valid sets the standard deviation varies by a factor of two ([Figure 25.4.2](#fig-cau-adjustment-sets)).
:::

\begin{center}
\begin{tikzpicture}[>=Stealth, line cap=round]
\definecolor{ink}{HTML}{1B1F27}\definecolor{accent}{HTML}{2F5F96}\definecolor{draft}{HTML}{A45C25}\definecolor{third}{HTML}{14705F}
\tikzset{var/.style={draw=ink, circle, minimum size=8mm, inner sep=0pt, font=\small}}
\node[var] (W1) at (0,2.2) {$W_1$};
\node[var] (W2) at (4.4,2.2) {$W_2$};
\node[var] (W3) at (2.2,1.25) {$W_3$};
\node[var] (I) at (-1.6,0.7) {$I$};
\node[var] (P) at (6.0,0.7) {$P$};
\node[var] (X) at (0,0) {$X$};
\node[var] (M) at (2.2,-1.0) {$M$};
\node[var] (Y) at (4.4,0) {$Y$};
\draw[->, thick, draft] (W1) -- (W2);
\draw[->, thick, draft] (W1) -- (X);
\draw[->, thick, draft] (W2) -- (Y);
\draw[->, thick, draft] (W3) -- (X);
\draw[->, thick, draft] (W3) -- (Y);
\draw[->, thick, ink] (I) -- (X);
\draw[->, thick, ink] (P) -- (Y);
\draw[->, thick, accent] (X) -- (Y);
\draw[->, thick, accent] (X) -- (M);
\draw[->, thick, accent] (M) -- (Y);
\end{tikzpicture}
\end{center}

[**Figure 25.4.1.** The graph of @exm-cau-eight. Blue arrows form the directed paths from \( X \) to \( Y \); brown arrows form
the two back-door paths. \( I \) influences only \( X \), and \( P \) only \( Y \).]{#fig-cau-eight-dag}

```{.python .run #cell-backdoor-model}
import numpy as np

names = ["W1", "W3", "I", "P", "W2", "X", "M", "Y"]      # a topological order
k = {v: j for j, v in enumerate(names)}
B = np.zeros((8, 8))                                      # B[j, i]: coefficient of i in j's equation
for child, parent, coef in [("W2", "W1", 0.8), ("X", "W1", 0.7), ("X", "W3", 0.6), ("X", "I", 0.9),
                            ("M", "X", 0.8), ("Y", "X", 0.5), ("Y", "M", 0.6), ("Y", "W2", 0.7),
                            ("Y", "W3", 0.8), ("Y", "P", 1.0)]:
    B[k[child], k[parent]] = coef
A = np.linalg.inv(np.eye(8) - B)                          # V = A U
Sigma = A @ A.T                                           # unit disturbance variances
tau = A[k["Y"], k["X"]]                                   # total effect: sum over directed paths

def population_fit(Z):
    """Coefficient of X in L(Y | X, Z) and the large-sample sd factor sigma_{Y.XZ}/sigma_{X.Z}."""
    R = [k["X"]] + [k[z] for z in Z]
    coef = np.linalg.solve(Sigma[np.ix_(R, R)], Sigma[R, k["Y"]])
    res_y = Sigma[k["Y"], k["Y"]] - Sigma[k["Y"], R] @ coef
    Zi = [k[z] for z in Z]
    res_x = Sigma[k["X"], k["X"]] - (Sigma[k["X"], Zi] @ np.linalg.solve(Sigma[np.ix_(Zi, Zi)], Sigma[Zi, k["X"]])
                                     if Z else 0.0)
    return coef[0], np.sqrt(res_y / res_x)

sets = [[], ["W3"], ["W1", "W3"], ["W2", "W3"], ["W2", "W3", "P"], ["W1", "W3", "I"],
        ["W1", "W3", "M"], ["W1", "W2", "W3", "I", "P"]]
print(f"total effect of X on Y: {tau:.3f}")
for Z in sets:
    c, f = population_fit(Z)
    print(f"Z = {{{', '.join(Z)}}}".ljust(26) + f"coefficient {c:6.3f}   sd factor {f:5.3f}")
```

::: {when-format="html"}
![**Figure 25.4.2.** The coefficient of \( X \) given each candidate set in @exm-cau-eight: population value (dot), a band of
\( \pm1.96 \) large-sample standard deviations at \( n=500 \), and the mean over the simulated samples (cross). Blue sets satisfy the
back-door criterion; the dashed line is \( \tau=0.98 \).](adjustment_sets.svg){#fig-cau-adjustment-sets width=80%}
:::

::: {when-format="pdf"}
![The coefficient of \( X \) given each candidate set in @exm-cau-eight: population value (dot), a band of
\( \pm1.96 \) large-sample standard deviations at \( n=500 \), and the mean over the simulated samples (cross). Blue sets satisfy the
back-door criterion; the dashed line is \( \tau=0.98 \).](adjustment_sets.pdf){width=80%}
:::

## Which adjustment set is most precise

All valid sets estimate \( \tau \), so the choice among them is a question of variance. With normal disturbances, \( Y \) given
\( (X,\mathbf{Z}) \) follows a linear model with variance \( \sigma^2_{Y\cdot XZ} \), and given the regressors the coefficient of \( X \)
has variance
\( \sigma^2_{Y\cdot XZ}/\norm{(\I-\M_Z)\x}^2 \), where \( \M_Z \) projects onto the span of \( \bone \) and the columns of \( \mathbf{Z} \) (@thm-proj-fwl and @eq-proj-vif-preview). Divided by \( n \), the squared length in the denominator is the residual mean square
of the regression of \( X \) on \( \mathbf{Z} \), which converges to \( \sigma^2_{X\cdot Z} \) with probability one (@prp-proj-consistency). So \( n \) times the variance tends to
\[
v(\mathcal Z)=\frac{\sigma^2_{Y\cdot XZ}}{\sigma^2_{X\cdot Z}},
\]{#eq-cau-asymptotic-variance}

the square of the factor in the table. A good adjustment set explains much of \( Y \) and little of \( X \).

::: {#prp-cau-adjustment-precision}
[Adding a variable to an adjustment set]

Let \( \mathcal Z \) and \( \mathcal Z\cup\{W\} \) both satisfy the back-door criterion, with all covariance matrices positive
definite.

::: {.enumerate options="label=(\alph*)"}
1. If \( \sigma_{YW\cdot XZ}=0 \), then \( v(\mathcal Z\cup\{W\})\ge v(\mathcal Z) \).

2. If \( \sigma_{XW\cdot Z}=0 \), then \( v(\mathcal Z\cup\{W\})\le v(\mathcal Z) \).
:::

:::

::: {.proof}
Partial variances are the mean squared errors of best linear predictors, and projecting onto a larger space can only
reduce them (@thm-rv-blp(c)). (a) By @prp-cor-partial-coefficient applied with \( (X,\mathbf{Z}) \) as conditioning variables,
\( \sigma_{YW\cdot XZ}=0 \) means that \( W \) has coefficient zero in \( L(Y\mid X,\mathbf{Z},W) \), so this predictor equals
\( L(Y\mid X,\mathbf{Z}) \) and the numerator of @eq-cau-asymptotic-variance is unchanged; the denominator can only decrease.
(b) Likewise \( \sigma_{XW\cdot Z}=0 \) means \( L(X\mid\mathbf{Z},W)=L(X\mid\mathbf{Z}) \), so the denominator is unchanged, while the
numerator can only decrease.
:::

In @exm-cau-eight, a direct computation from \( \Cov(\mathbf{V})=\A\A\T \) gives \( \sigma_{YI\cdot XW_1W_3}=0 \) and
\( \sigma_{XP\cdot W_2W_3}=0 \): the pure cause \( I \) of \( X \) satisfies (a) and the pure cause \( P \) of \( Y \) satisfies (b).
Accordingly, adding \( I \) to \( \{W_1,W_3\} \) raises the factor from \( 1.25 \) to
\( 1.69 \); adding \( P \) to \( \{W_2,W_3\} \) lowers it from \( 1.06 \)
to \( 0.80 \). Replacing \( W_1 \) by \( W_2 \), closer to \( Y \) on the back-door path, also helps. Henckel, Perković and
Maathuis (2022) turn these observations into a graphical rule for the valid set of smallest asymptotic variance in linear
models: roughly, the parents of \( Y \) and of the mediators, other than the mediators and \( X \).

If every adjustment set contains an unrecorded variable, as in @exm-cau-visits, other structure must be used: an instrumental variable (@thm-eiv-iv in
[Chapter 24](../ch24-errors-in-variables/index.html)), or the front-door criterion of Pearl (1995), which uses an observed
mediator; or a sensitivity analysis asks how strong an unobserved confounder would have to be, using formulas such as @eq-cau-backdoor-bias.

## Exercises

### A. Check your understanding

::: {#exr-cau-three-sets}
[A1]

In @exm-cau-eight, decide which of \( \{W_1,W_2\} \), \( \{W_3,P\} \) and \( \{W_2,W_3,M\} \) satisfy the back-door
criterion, and explain the coefficients of \( X \) that the listing gives for them: \( 1.201 \), \( 1.150 \) and
\( 0.500 \).
:::

::: {.solution}
None does. \( \{W_1,W_2\} \) leaves \( X\leftarrow W_3\to Y \) open, and \( \{W_3,P\} \) leaves
\( X\leftarrow W_1\to W_2\to Y \) open, so both coefficients exceed \( \tau=0.98 \). The set \( \{W_2,W_3,M\} \) blocks both
back-door paths but contains the mediator \( M \), violating (i); holding \( M \) fixed leaves only the direct arrow, with
coefficient \( 0.5 \).
:::

### B. Practice

::: {#exr-cau-stratified}
[B1]

Let \( X\in\{0,1\} \) and let \( Z \) take values \( 1,\dots,K \). In a sample, stratum \( z \) has \( n_z \) units, a fraction
\( \hat{p}_z\in(0,1) \) of them with \( X=1 \), and difference of mean responses \( \hat{\delta}_z=\bar{y}_{1z}-\bar{y}_{0z} \). Show that the sample
version of @eq-cau-adjustment for the effect \( \E\bigl(Y\mid\operatorname{do}(X=1)\bigr)-\E\bigl(Y\mid\operatorname{do}(X=0)\bigr) \) is
\( \sum_z(n_z/n)\hat{\delta}_z \), while the coefficient of \( X \) in the least squares regression of \( Y \) on \( X \) and the \( K \) stratum
indicators is
\[
\hat{\beta}_X=\frac{\sum_zn_z\hat{p}_z(1-\hat{p}_z)\,\hat{\delta}_z}{\sum_zn_z\hat{p}_z(1-\hat{p}_z)}.
\]
When do they agree?
:::

::: {.solution}
Replacing \( \E(Y\mid X=x,Z=z) \) by \( \bar{y}_{xz} \) and the distribution of \( Z \) by the sample fractions \( n_z/n \) gives the first
expression. For the second, apply @thm-proj-fwl with \( \X_2 \) the stratum indicators. Projecting onto their span replaces each
value by its stratum mean (@exm-proj-oneway-M), so \( \tilde{x}_i=x_i-\hat{p}_{z(i)} \). Then
\( \hat{\beta}_X=\sum_i\tilde{x}_iy_i/\sum_i\tilde{x}_i^2 \). Within stratum \( z \),
\( \sum_i\tilde{x}_i^2=n_z\hat{p}_z(1-\hat{p}_z) \) and
\( \sum_i\tilde{x}_iy_i=n_z\hat{p}_z(1-\hat{p}_z)(\bar{y}_{1z}-\bar{y}_{0z}) \), since the treated contribute
\( (1-\hat{p}_z)n_z\hat{p}_z\bar{y}_{1z} \) and the untreated \( -\hat{p}_zn_z(1-\hat{p}_z)\bar{y}_{0z} \). The two agree when the
\( \hat{\delta}_z \) are all equal or when \( \hat{p}_z \) is the same in every stratum. Otherwise the regression weights strata by the
variance of \( X \) within them, not by their size, and estimates a weighted average effect that is not the average effect
(Angrist 1998).
:::

::: {#exr-cau-table-two}
[B2]

In @exm-cau-eight, show that the coefficient of \( W_3 \) in \( L(Y\mid X,W_1,W_3) \) is \( 0.80 \), the coefficient of the
arrow \( W_3\to Y \), while the total effect of \( W_3 \) on \( Y \) is \( 1.388 \). Which set of variables would one adjust for to
estimate the effect of \( W_3 \)?
:::

::: {.solution}
Write \( Y=0.5X+0.6(0.8X+U_M)+0.7(0.8W_1+U_{W_2})+0.8W_3+P+U_Y=0.98X+0.56W_1+0.8W_3+E \), where
\( E=0.6U_M+0.7U_{W_2}+P+U_Y \) is uncorrelated with \( X \), \( W_1 \) and \( W_3 \) (none of these depends on the disturbances in \( E \)). So
\( L(Y\mid X,W_1,W_3)=0.98X+0.56W_1+0.8W_3 \). The total effect of \( W_3 \) is the sum over its directed paths to \( Y \):
\( 0.8 \) directly plus \( 0.6\times0.98 \) through \( X \), which is \( 1.388 \). The coefficient in the regression misses the part
that goes through \( X \), because \( X \) is held fixed. \( W_3 \) has no parents and no back-door paths, so the empty set is an
adjustment set for its effect: the simple regression of \( Y \) on \( W_3 \) estimates \( 1.388 \).
:::

::: {#exr-cau-randomized-precision}
[B3]

Suppose \( X \) is randomized, so it has no parents, and \( P \) is a pre-treatment cause of \( Y \). Show that both \( \emptyset \) and
\( \{P\} \) are adjustment sets and that \( v(\{P\})=(1-\rho^2_{YP\cdot X})\,v(\emptyset) \), where \( \rho_{YP\cdot X} \) is the partial
correlation of \( Y \) and \( P \) given \( X \). Compare with @prp-dsn-precision.
:::

### C. Going deeper

::: {#exr-cau-positivity}
[C1]

Positivity is not a technicality. Let \( Z \) be \( 1 \) or \( 0 \) with probability \( \tfrac12 \) each; let
\( X \) be \( 1 \) or \( 0 \) with probability \( \tfrac12 \) each when \( Z=0 \), and \( X=0 \) whenever
\( Z=1 \); and let
\[
Y=\tau X+\lambda Z+\psi XZ+U_Y,
\]
with \( U_Y \) of mean zero and independent of \( (X,Z) \). The graph is \( Z\to X \), \( Z\to Y \),
\( X\to Y \), so \( \{Z\} \) satisfies the back-door criterion.

::: {.enumerate options="label=(\alph*)"}
1. Show that the average effect
   \( \E\bigl(Y\mid\operatorname{do}(X=1)\bigr)-\E\bigl(Y\mid\operatorname{do}(X=0)\bigr) \) is
   \( \tau+\psi/2 \).

2. Show that the right-hand side of @eq-cau-adjustment is not defined, because the positivity condition of
   @thm-cau-backdoor(a) fails at \( z=1 \).

3. Show that the joint distribution of \( (X,Z,Y) \) does not involve \( \psi \) at all, and that the least
   squares regression of \( Y \) on \( 1 \), \( X \) and \( Z \) nevertheless returns the coefficient
   \( \tau \) exactly, in the population and in every sample in which the three covariate patterns occur.
:::

What, then, is the regression's answer an answer to?
:::

::: {.solution}
(a) Under \( \operatorname{do}(X=x) \) the assignment for \( Z \) is unchanged, so
\( \E\bigl(Y\mid\operatorname{do}(X=x)\bigr)=\tau x+\lambda\E Z+\psi x\E Z=\tau x+(\lambda+\psi x)/2 \), and
the difference at \( x=1 \) and \( x=0 \) is \( \tau+\psi/2 \).

(b) \( \Pr(X=1\mid Z=1)=0 \), so \( \E(Y\mid X=1,Z=1) \) is a conditional expectation on a null event: the
distribution of the data says nothing about it, and the outer expectation in @eq-cau-adjustment averages a
quantity that is undefined on half of the population.

(c) The event \( \{X=1,Z=1\} \) has probability zero, so \( XZ=0 \) with probability one and
\( Y=\tau X+\lambda Z+U_Y \) almost surely: no feature of the distribution moves when \( \psi \) changes, and
\( \psi \) is not identified. The covariate vector \( (1,X,Z) \) takes only the three values
\( (1,0,0) \), \( (1,1,0) \) and \( (1,0,1) \), which are linearly independent, so the linear model is
saturated: the best linear predictor reproduces \( \E(Y\mid X,Z) \) at each of the three patterns, giving
\( \beta_0=0 \), \( \beta_0+\beta_1=\tau \) and \( \beta_0+\beta_2=\lambda \). Hence \( \beta_1=\tau \),
and in a sample the normal equations have the same unique solution in terms of the three cell means.

So the regression reports \( \tau \), which is the effect of \( X \) among the units with \( Z=0 \), and
equals the average effect only under the assumption \( \psi=0 \) — an assumption the data cannot examine,
because the stratum in which \( X=1 \) and \( Z=1 \) is empty. The number comes from the functional form, by
extrapolation, and not from the observations. A linear model is at its most dangerous where the design is
thin, since it never declines to answer.
:::

::: {#exr-cau-front-door}
[C2]

When every adjustment set contains an unrecorded variable, an observed *mediator* can still identify the
effect. Let \( U \) be unrecorded and put
\[
X=aU+U_X,\qquad M=\delta X+U_M,\qquad Y=\gamma M+bU+U_Y,
\]
with \( U,U_X,U_M,U_Y \) independent of mean zero and variances
\( \sigma_U^2,\sigma_X^2,\sigma_M^2,\sigma_Y^2 \), all positive, and \( a\ne0 \), \( b\ne0 \).

::: {.enumerate options="label=(\alph*)"}
1. Show that neither \( \emptyset \) nor \( \{M\} \) satisfies the back-door criterion relative to
   \( (X,Y) \), so @thm-cau-backdoor identifies nothing from the recorded variables \( (X,M,Y) \).

2. Show that the slope of \( L(M\mid X) \) is exactly \( \delta \).

3. Show that the coefficient of \( M \) in \( L(Y\mid X,M) \) is exactly \( \gamma \), while the coefficient
   of \( X \) there is \( b\kappa \) with \( \kappa=a\sigma_U^2/(a^2\sigma_U^2+\sigma_X^2) \), and so is not
   the direct effect \( 0 \).

4. Conclude that the product of the two estimable slopes is the total effect \( \tau=\gamma\delta \). This is
   the linear case of Pearl's front-door criterion.

5. Show that the construction is fragile: if \( M=\delta X+gU+U_M \) with \( g\ne0 \), the two slopes become
   \( \delta+g\kappa \) and \( \gamma+b\lambda \) with
   \( \lambda=g\sigma_U^2/(g^2\sigma_U^2+\sigma_M^2) \), whose product is not \( \gamma\delta \).
:::

:::

::: {.solution}
(a) The path \( X\leftarrow U\to Y \) is a back-door path with \( U \) a non-collider, so the empty set fails
(ii); and \( M \) is a descendant of \( X \), so \( \{M\} \) fails (i). Since \( U \) is unrecorded, there is
no candidate set left.

(b) \( \Cov(M,X)=\delta\Var(X)+\Cov(U_M,X)=\delta\Var(X) \), because \( U_M \) is independent of \( X \), so
@thm-proj-blp gives the slope \( \Cov(M,X)/\Var(X)=\delta \).

(c) Since \( U_Y \) is uncorrelated with \( (X,M) \), \( L(Y\mid X,M)=\gamma M+b\,L(U\mid X,M) \). Now
\( M=\delta X+U_M \), so the span of \( 1,X,M \) is the span of \( 1,X,U_M \), and \( U_M \) is uncorrelated
with \( X \); by the orthogonality of the two pieces,
\( L(U\mid X,M)=L(U\mid X)+L(U\mid U_M)=\kappa X+0 \), with
\( \kappa=\Cov(U,X)/\Var(X)=a\sigma_U^2/(a^2\sigma_U^2+\sigma_X^2) \). So
\( L(Y\mid X,M)=\gamma M+b\kappa X \): the coefficient of \( M \) is exactly \( \gamma \), and the
coefficient of \( X \) is \( b\kappa\ne0 \) although the arrow \( X\to Y \) is absent. Conditioning on the
mediator has left the confounding in the coefficient of \( X \) and removed it from the coefficient of
\( M \) — because every back-door path from \( M \) to \( Y \) runs through \( X \), which is held fixed.

(d) The total effect of \( X \) on \( Y \) is the product over the one directed path, \( \tau=\gamma\delta \),
and by (b) and (c) both factors are coefficients of regressions among recorded variables, so each is
estimated consistently by least squares (@prp-proj-consistency).

(e) Replace \( U_M \) by \( E=gU+U_M \), which is still uncorrelated with \( X \), and repeat the two
computations. For (b), \( \Cov(M,X)=\delta\Var(X)+g\,a\sigma_U^2 \), so the slope is \( \delta+g\kappa \).
For (c), \( L(U\mid X,M)=L(U\mid X)+L(U\mid E)=\kappa X+\lambda E \) with
\( \lambda=g\sigma_U^2/(g^2\sigma_U^2+\sigma_M^2) \), and \( E=M-\delta X \), so the coefficient of \( M \)
is \( \gamma+b\lambda \). The product \( (\delta+g\kappa)(\gamma+b\lambda) \) differs from \( \gamma\delta \)
by terms of first order in \( g \). The front-door criterion buys identification from an assumption — that
nothing unrecorded acts on the mediator — that is no more checkable than the one it replaces.
:::
