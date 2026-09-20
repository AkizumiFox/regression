# Deviance and asymptotic inference

Parts II and III rested on exact distribution theory: \( \hbeta \) was normal because \( \Y \)
was, and sums of squares were chi-squared because they were quadratic forms in a normal
vector (@thm-qf-chisq). None of that survives. This section states what replaces it in the
limit, under which conditions, and — just as important — where the approximations fail.

## Asymptotic normality

The conditions below are stated for the canonical link, where the argument can be given in
full. Throughout \( \bbeta^0 \) is the true parameter, \( \theta_i^0=\x_{(i)}\T\bbeta^0 \), and
\[
\A_n=\boldsymbol{\mathcal I}_n(\bbeta^0)=\frac1\phi\,\X\T\diag\{w_i\,b''(\theta_i^0)\}\X
\]
is the expected information from @eq-glm-information. Because the canonical-link Hessian
does not involve \( \y \) (@thm-glm-score(d)), \( \A_n \) is a fixed matrix, not a random one.
Define the scalars
\[
q_{ni}=\frac{w_i\,b''(\theta_i^0)}{\phi}\;\x_{(i)}\T\A_n^{-1}\x_{(i)},\qquad i=1,\dots,n .
\]{#eq-glm-leverage}

These are the exact analogue of leverages: nonnegative, summing
to \( \tr(\A_n^{-1}\A_n)=p \), so they average \( p/n \) as
in @prp-proj-leverage(b). The conditions are

- **(G1)** the responses are independent, the model is that of @def-glm-model with the
  canonical link, \( \rank(\X)=p \) with \( p \) fixed, and \( \phi \) and the \( w_i \) are known;
- **(G2)** \( \bbeta^0 \) is interior to \( \mathcal B \), the natural parameters \( \theta_i^0 \)
  lie in a fixed compact subset \( K \) of the interior of \( \Theta \), and
  \( 0<w_{\min}\le w_i\le w_{\max}<\infty \);
- **(G3)** \( \max_{i\le n}q_{ni}\to0 \): no single observation carries a fixed share of the
  information;
- **(G4)** \( \lambda_{\min}(\A_n)\to\infty \): information accumulates in every direction.

::: {#thm-glm-asymptotics}
[Consistency and asymptotic normality of the estimate]

Under (G1)–(G4), with probability tending to one the likelihood equations have a solution
\( \hbeta_n \); it is the unique maximizer of \( \ell \), it is consistent for \( \bbeta^0 \), and
\[
\A_n^{1/2}\bigl(\hbeta_n-\bbeta^0\bigr)\ \xrightarrow{d}\ \Normal_p(\bzero,\I).
\]{#eq-glm-asymptotic-normality}

The same holds with \( \A_n \) replaced by \( \hat{\A}_n=\boldsymbol{\mathcal I}_n(\hbeta_n) \).
Consequently, for a fixed \( \blambda\ne\bzero \),
\[
\frac{\blambda\T\hbeta_n-\blambda\T\bbeta^0}{\bigl(\blambda\T\hat{\A}_n^{-1}\blambda\bigr)^{1/2}}
\ \xrightarrow{d}\ \Normal(0,1),
\]
and \( \blambda\T\hbeta_n\pm z_{1-\alpha/2}(\blambda\T\hat{\A}_n^{-1}\blambda)^{1/2} \) is an
interval with asymptotic coverage \( 1-\alpha \).
:::

One result is imported, not proved: the standard corollary of **the convexity lemma** of
Pollard (1991, section 6), that pointwise convergence of concave random functions is
automatically uniform on compact sets. Let \( \Lambda_n \) be random concave functions on
\( \Real^p \), allowed the value \( -\infty \) outside a convex set, of the form
\( \Lambda_n(\mathbf{s})=\mathbf{s}\T\bz_n-\tfrac12\norm{\mathbf{s}}^2+r_n(\mathbf{s}) \) with
\( \bz_n=O_p(1) \) and \( \sup_{\norm{\mathbf{s}}\le M}\lvert r_n(\mathbf{s})\rvert\to0 \) in
probability for every \( M \). Then a maximizer \( \hat{\mathbf{s}}_n \) exists with probability
tending to one and \( \hat{\mathbf{s}}_n-\bz_n\to\bzero \) in probability. Concavity does
the work that compactness conditions do in the usual theory, and @thm-glm-concave supplies
it free.

::: {.proof}
*Step 1: the score is asymptotically normal.* By @eq-glm-score-canonical,
\( \bU_n=\phi^{-1}\sum_iw_i(Y_i-\mu_i^0)\x_{(i)} \) is a sum of independent mean-zero vectors
with \( \Cov(\bU_n)=\A_n \). Put \( \bz_n=\A_n^{-1/2}\bU_n \). For a fixed unit vector
\( \bv \) write \( \bv\T\bz_n=\sum_ic_{ni}\xi_i \) with \( \xi_i=Y_i-\mu_i^0 \) and the
weights \( c_{ni}=\phi^{-1}w_i\,\bv\T\A_n^{-1/2}\x_{(i)} \). This is the weighted sum
of @lem-het-clt, and its case (b) applies. Moments: by (G2) the cumulants of \( Y_i \) are
continuous functions of \( \theta_i^0 \) on a compact set (@prp-glm-moments(c)) and the
\( w_i \) are bounded away from \( 0 \) and \( \infty \), so \( \E\xi_i^4 \) is bounded — case (b)
with \( \delta=2 \) — and \( \tau_{\min}^2\le\Var\xi_i=\phi b''(\theta_i^0)/w_i\le\tau_{\max}^2 \)
for constants free of \( i \) and \( n \). Weights: by construction
\( \sum_ic_{ni}^2\Var\xi_i=\bv\T\A_n^{-1/2}\A_n\A_n^{-1/2}\bv=1 \), so the sum is
already normalized and \( \norm{\mathbf{c}_n}^2\ge\tau_{\max}^{-2} \); while
\( c_{ni}^2\le\phi^{-2}w_{\max}^2\norm{\A_n^{-1/2}\x_{(i)}}^2 \) and
\( \norm{\A_n^{-1/2}\x_{(i)}}^2=\x_{(i)}\T\A_n^{-1}\x_{(i)}=\phi\,q_{ni}/\{w_ib''(\theta_i^0)\}\le C'q_{ni} \)
by @eq-glm-leverage and (G2) again. Hence
\( \max_ic_{ni}^2/\norm{\mathbf{c}_n}^2\le C''\max_iq_{ni}\to0 \) by (G3), which
is @eq-het-no-dominant-weight. So \( \bv\T\bz_n\to\Normal(0,1) \) for every unit
\( \bv \), and the Cramér–Wold device gives \( \bz_n\to\Normal_p(\bzero,\I) \).

*Step 2: the log-likelihood is asymptotically quadratic.* Put
\( \Lambda_n(\mathbf{s})=\ell(\bbeta^0+\A_n^{-1/2}\mathbf{s})-\ell(\bbeta^0) \). Taylor's theorem
with integral remainder, applied to the smooth function \( \ell \), gives
\[
\Lambda_n(\mathbf{s})=\mathbf{s}\T\bz_n
-\int_0^1(1-t)\,\mathbf{s}\T\A_n^{-1/2}\boldsymbol{\mathcal J}_n(\bbeta_t)\A_n^{-1/2}\mathbf{s}\,dt,
\]
where \( \bbeta_t=\bbeta^0+t\A_n^{-1/2}\mathbf{s} \) and
\( \boldsymbol{\mathcal J}_n(\bbeta)=\phi^{-1}\X\T\diag\{w_ib''(\x_{(i)}\T\bbeta)\}\X \) is the
*nonrandom* observed information. Fix \( M>0 \) and let \( \norm{\mathbf{s}}\le M \). Then
\( \lvert\x_{(i)}\T\A_n^{-1/2}\mathbf{s}\rvert\le M\norm{\A_n^{-1/2}\x_{(i)}}\le M(C'q_{ni})^{1/2} \),
which by (G3) tends to \( 0 \) uniformly in \( i \). So for large \( n \) every perturbed natural
parameter \( \theta_i^0+t\,\x_{(i)}\T\A_n^{-1/2}\mathbf{s} \) lies in a fixed compact
neighbourhood of \( K \) inside \( \Theta \), on which \( b'' \) is uniformly continuous and
bounded away from \( 0 \) and \( \infty \); hence
\( b''(\theta_i^0+\delta_i)=b''(\theta_i^0)\{1+\epsilon_{ni}\} \) with
\( \sup_i\lvert\epsilon_{ni}\rvert\to0 \) uniformly in \( t\in[0,1] \) and
\( \norm{\mathbf{s}}\le M \), and
\( \A_n^{-1/2}\boldsymbol{\mathcal J}_n(\cdot)\A_n^{-1/2}=\I+\R_n \) with \( \norm{\R_n}\to0 \), and
\[
\Lambda_n(\mathbf{s})=\mathbf{s}\T\bz_n-\tfrac12\norm{\mathbf{s}}^2+r_n(\mathbf{s}),
\qquad \sup_{\norm{\mathbf{s}}\le M}\lvert r_n(\mathbf{s})\rvert\to0 .
\]{#eq-glm-quadratic-expansion}

*Step 3: the maximizer.* By @thm-glm-concave each \( \Lambda_n \) is concave. It is defined
only on the convex set \( \A_n^{1/2}(\mathcal B-\bbeta^0) \), a half-space for the canonical
links of the gamma and inverse Gaussian; extend it by \( -\infty \) outside, which preserves
concavity, and note that Step 2 put every ball \( \norm{\mathbf{s}}\le M \) inside that set for
large \( n \), so the expansion holds where the lemma needs it.
Then @eq-glm-quadratic-expansion is exactly the form the lemma requires, with
\( \bz_n=O_p(1) \) by Step 1 and a nonrandom remainder bound, so a maximizer
\( \hat{\mathbf{s}}_n \) exists with probability tending to one and
\( \hat{\mathbf{s}}_n-\bz_n\to\bzero \) in probability. Since
\( \hat{\mathbf{s}}_n=\A_n^{1/2}(\hbeta_n-\bbeta^0) \), Slutsky's theorem and Step 1
give @eq-glm-asymptotic-normality, and @thm-glm-concave makes the maximizer unique.

*Consistency.* \( \norm{\hbeta_n-\bbeta^0}\le\lambda_{\min}(\A_n)^{-1/2}\norm{\hat{\mathbf{s}}_n} \),
and \( \norm{\hat{\mathbf{s}}_n}=O_p(1) \), so (G4) gives \( \hbeta_n\to\bbeta^0 \) in probability.

*Plug-in.* Consistency and the uniform continuity of \( b'' \) used in Step 2 give
\[
\A_n^{-1/2}\hat{\A}_n\A_n^{-1/2}\to\I
\]
in probability, whence \( \hat{\A}_n^{1/2}(\hbeta_n-\bbeta^0)\to\Normal_p(\bzero,\I) \). The
last two displays follow by taking the \( \blambda \) component.
:::

::: {.remark}
[What the conditions mean, and the non-canonical case]

Condition (G3) is the honest content of "large sample". It is the generalized-linear-model
version of the requirement \( \max_ih_{ii}\to0 \) that
[Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) needed, and Step 1
used it in exactly that role, to supply @eq-het-no-dominant-weight. It fails when it should:
a group with a single observation, a regressor that indicates one case, a binary design with
a nearly empty cell. The canonical link is used twice, for the concavity and for the
nonrandom Hessian, and the next proposition records what survives without it.
:::

::: {.warning}
Nothing here is a finite-sample statement: \( \hbeta \) is biased, by \( O(n^{-1}) \) against a
standard error of \( O(n^{-1/2}) \) — asymptotically negligible, not zero. Nor does \( p \)
fixed cover many parameters: a model with a parameter per group and bounded group sizes
violates (G3), and @exm-opt-neyman-scott showed what goes wrong.
:::

**Intervals for a fitted mean.** For fixed \( \x_0 \) the theorem with \( \blambda=\x_0 \)
gives \( \hat\eta_0\pm z_{1-\alpha/2}(\x_0\T\hat{\A}_n^{-1}\x_0)^{1/2} \) for
\( \eta_0=\x_0\T\bbeta^0 \), the counterpart of @thm-ci-estimable-interval with a normal
quantile for a \( t \) one; applying the monotone \( h \) to both endpoints makes it an
interval for \( \mu_0=h(\eta_0) \) of the same asymptotic coverage. The *interval*
transforms, not the endpoint formula, so it is not symmetric about \( \hat\mu_0 \) — a
feature, since a log link then keeps it positive and a logit link inside \( (0,1) \). The
delta-method alternative
\( \hat\mu_0\pm z_{1-\alpha/2}\,h'(\hat\eta_0)(\x_0\T\hat{\A}_n^{-1}\x_0)^{1/2} \) agrees
to first order and is the worse wherever \( h \) curves sharply. Neither is a prediction
interval: @thm-ci-prediction-interval has no clean analogue here.

::: {#prp-glm-noncanonical}
[Non-canonical links]

Let the link be any twice continuously differentiable \( g \) with \( g'\ne0 \) on
\( \mathcal M \), and put \( \A_n=\phi^{-1}\X\T\W\X \) with \( \W \) the working weights
of @thm-glm-score evaluated at \( \bbeta^0 \). Then every conclusion
of @thm-glm-asymptotics — existence with probability tending to one, consistency, the
limit @eq-glm-asymptotic-normality, the plug-in version and the Wald intervals — continues to
hold, as do @prp-glm-three-tests and @thm-glm-deviance(a),(b) below, which use the theorem
only through the quadratic expansion @eq-glm-quadratic-expansion. The conditions are those
of Fahrmeir and Kaufmann (1985): the divergence condition
\( \lambda_{\min}(\A_n)\to\infty \), and a continuity condition asking that
\( \A_n(\bbeta^0)^{-1/2}\A_n(\bbeta)\A_n(\bbeta^0)^{-1/2}\to\I \) uniformly over
\( \bbeta \) in neighbourhoods of \( \bbeta^0 \) that shrink in the information metric. The
two together play the part of (G2)–(G4).
:::

::: {.remark}
This proposition is quoted, not proved: without concavity and a nonrandom Hessian the
argument above collapses, and the replacement is the content of Fahrmeir and
Kaufmann (1985). It is recorded because most models in the chapters that follow use a
non-canonical link somewhere — the complementary log–log and probit links
of [Chapter 35](../ch35-binary-responses/index.html), the gamma/log fit
of @exm-glm-engel, the negative binomial with a log link
in [Chapter 37](../ch37-counts/index.html).
:::

## Wald, score and likelihood ratio

Consider a reduced model \( \C(\X_0)\subseteq\C(\X) \) with \( \rank(\X_0)=p-q \), and let
\( \tilde{\bbeta} \) be the maximum likelihood estimate under it, written as a vector in
\( \Real^p \). Since \( \X \) has full column rank, \( \{\bbeta:\X\bbeta\in\C(\X_0)\} \) is a
subspace of dimension \( p-q \), so it is \( \{\bbeta:\bLambda\T\bbeta=\bzero\} \) for some
\( p\times q \) matrix \( \bLambda \) of rank \( q \), in the notation
of @def-glh-testable.

::: {#prp-glm-three-tests}
[The three statistics]

Define
\[
\begin{aligned}
\text{Wald}&=(\bLambda\T\hbeta)\T\bigl(\bLambda\T\hat{\A}_n^{-1}\bLambda\bigr)^{-1}\bLambda\T\hbeta,\\
\text{Score}&=\bU(\tilde{\bbeta})\T\boldsymbol{\mathcal I}_n(\tilde{\bbeta})^{-1}\bU(\tilde{\bbeta}),\\
\text{LR}&=2\bigl\{\ell(\hbeta)-\ell(\tilde{\bbeta})\bigr\}.
\end{aligned}
\]
Under (G1)–(G4) and the reduced model, each of the three converges in distribution to
\( \chi^2(q) \), and any two of them differ by \( o_p(1) \). For a single coefficient
(\( q=1 \), \( \bLambda=\mathbf{e}_j \)) the Wald statistic is
\( \hat\beta_j^2/[\hat{\A}_n^{-1}]_{jj} \), the square of the "estimate over standard error"
ratio that every program prints.
:::

::: {.proof}
Work in the coordinates \( \mathbf{s}=\A_n^{1/2}(\bbeta-\bbeta^0) \) of the proof
of @thm-glm-asymptotics, where \( \bbeta^0 \) is the true parameter, which lies in the reduced
model, so \( \bLambda\T\bbeta^0=\bzero \). The reduced model is the linear subspace
\( \mathcal S_n=\A_n^{1/2}\Null(\bLambda\T)=\Null(\B_n) \) of \( \Real^p \), where
\( \B_n=\bLambda\T\A_n^{-1/2} \) has rank \( q \); let \( \bP_n \) be the orthogonal projection
onto it, so that \( \I-\bP_n=\B_n\T(\B_n\B_n\T)^{-1}\B_n \).

For the likelihood ratio, @eq-glm-quadratic-expansion gives
\( \Lambda_n(\mathbf{s})=\mathbf{s}\T\bz_n-\tfrac12\norm{\mathbf{s}}^2+r_n(\mathbf{s}) \) uniformly on
compacts. Maximizing \( \mathbf{s}\T\bz_n-\tfrac12\norm{\mathbf{s}}^2 \) over \( \Real^p \) gives
\( \tfrac12\norm{\bz_n}^2 \) at \( \mathbf{s}=\bz_n \); over \( \mathcal S_n \) it gives
\( \tfrac12\norm{\bP_n\bz_n}^2 \) at \( \mathbf{s}=\bP_n\bz_n \). Both maximizers are \( O_p(1) \), so
the remainder is negligible at both, and
\[
\text{LR}=2\bigl\{\max_{\Real^p}\Lambda_n-\max_{\mathcal S_n}\Lambda_n\bigr\}
=\norm{(\I-\bP_n)\bz_n}^2+o_p(1).
\]
If \( \bz\sim\Normal_p(\bzero,\I) \) then \( \norm{(\I-\bP)\bz}^2\sim\chi^2(q) \) for *every*
projection \( \bP \) of rank \( p-q \) (@thm-qf-chisq). Projections of a fixed rank form a compact
set, so along any subsequence \( \bP_n\to\bP \), and the continuous mapping theorem gives
\( \norm{(\I-\bP_n)\bz_n}^2\to\chi^2(q) \); the limit is the same for every subsequence, so the
whole sequence converges.

For the Wald statistic, write \( \hat{\mathbf{s}}_n=\A_n^{1/2}(\hbeta-\bbeta^0) \). Then
\( \bLambda\T\hbeta=\bLambda\T(\hbeta-\bbeta^0)=\B_n\hat{\mathbf{s}}_n \) and
\( \bLambda\T\A_n^{-1}\bLambda=\B_n\B_n\T \), so with \( \A_n \) in place of \( \hat{\A}_n \)
the statistic is exactly \( \norm{(\I-\bP_n)\hat{\mathbf{s}}_n}^2 \). The proof
of @thm-glm-asymptotics gives \( \hat{\mathbf{s}}_n=\bz_n+o_p(1) \), and
\( \A_n^{-1/2}\hat{\A}_n\A_n^{-1/2}\to\I \) in probability makes the substitution of
\( \hat{\A}_n \) cost \( o_p(1) \); hence \( \text{Wald}=\norm{(\I-\bP_n)\bz_n}^2+o_p(1) \). The
same argument inside the reduced model gives
\( \A_n^{1/2}(\tilde{\bbeta}-\bbeta^0)=\bP_n\bz_n+o_p(1) \), which is used next.

For the score statistic, differentiating @eq-glm-quadratic-expansion — legitimate because
pointwise convergence of differentiable concave functions forces convergence of their
gradients (Rockafellar 1970, section 25) — gives
\( \A_n^{-1/2}\bU(\bbeta^0+\A_n^{-1/2}\mathbf{s})=\bz_n-\mathbf{s}+o_p(1) \) uniformly on compacts;
at \( \mathbf{s}=\bP_n\bz_n+o_p(1) \) this is \( (\I-\bP_n)\bz_n+o_p(1) \).
:::

In the normal linear model the three are *exact* monotone functions of one another and of
the \( F \) statistic, with \( \text{Wald}>\text{LR}>\text{Score} \) in every
sample (@prp-glh-trinity). Outside it that ordering is not guaranteed, and the three agree
only to first order.

**The Wald statistic is not invariant.** Testing \( \beta_j=0 \) and testing
\( e^{\beta_j}=1 \) are the same hypothesis with different Wald statistics, because the
delta-method standard error of \( e^{\hat\beta_j} \) is not the transform of that of
\( \hat\beta_j \); the likelihood ratio statistic, which compares maximized likelihoods, is
invariant under every reparameterization.

**The Wald statistic can move the wrong way.** This failure has a name.

::: {#exm-glm-hauck-donner}
[The Hauck–Donner effect]

Two groups of \( 40 \) trials, with \( 20 \) successes in group A and \( s_B \) in group B, fitted
by a logistic model with an intercept and a group indicator. The coefficient is the log odds
ratio \( \hat\beta_1=\log\{s_B/(40-s_B)\} \), with standard error
\( \{1/20+1/20+1/s_B+1/(40-s_B)\}^{1/2} \). As \( s_B \) rises from \( 21 \) to \( 39 \) the
evidence against equality strengthens monotonically and the likelihood ratio statistic
records it, climbing to \( 27.301 \). The Wald statistic does not: it rises
to \( 13.842 \) at \( s_B=38 \) and then *falls* to
\( 11.924 \) at \( s_B=39 \), though the fitted log odds ratio has grown
to \( 3.664 \). The standard error grows faster than the coefficient.

At complete separation, \( s_B=40 \), coefficient and standard error are both infinite and
the Wald statistic is \( 0 \): the test accepts the null precisely where the data support it
most strongly
([Figure 34.5.1](05-deviance-and-inference.html#fig-glm-hauck-donner)). The effect occurs
whenever the log-likelihood is far from quadratic over the range where the estimate lies.
:::

::: {when-format="html"}
![**Figure 34.5.1.** The Hauck–Donner effect. Two groups of 40 trials, with 20 successes in
group A. As the successes in group B increase, the likelihood ratio statistic rises
monotonically, while the Wald statistic turns over and falls. The horizontal line is the
upper \( 5\% \) point of \( \chi^2(1) \).](hauck_donner.svg){#fig-glm-hauck-donner width=62%}
:::

::: {when-format="pdf"}
![The Hauck–Donner effect. Two groups of 40 trials, with 20 successes in
group A. As the successes in group B increase, the likelihood ratio statistic rises
monotonically, while the Wald statistic turns over and falls. The horizontal line is the
upper \( 5\% \) point of \( \chi^2(1) \).](hauck_donner.pdf){width=62%}
:::

```{.python .run #cell-inference-hauck-donner}
import numpy as np

m_hd, s_a = 40, 20                                  # group A: 20 successes out of 40
s_b = np.arange(21, 40)                             # group B: from 21 to 39 successes
odds_ratio = (s_b / (m_hd - s_b)) / (s_a / (m_hd - s_a))
beta1 = np.log(odds_ratio)                          # the fitted log odds ratio
se1 = np.sqrt(1 / s_a + 1 / (m_hd - s_a) + 1 / s_b + 1 / (m_hd - s_b))
wald = (beta1 / se1) ** 2

def binom_ll(s, m, p):
    return s * np.log(p) + (m - s) * np.log(1 - p)

p_pool = (s_a + s_b) / (2 * m_hd)
lr = 2 * (binom_ll(s_a, m_hd, s_a / m_hd) + binom_ll(s_b, m_hd, s_b / m_hd)
          - binom_ll(s_a, m_hd, p_pool) - binom_ll(s_b, m_hd, p_pool))

for k in (0, 10, 15, 17, 18):
    print(f"successes in B = {s_b[k]:2d}:  log odds ratio {beta1[k]:6.3f}"
          f"   Wald {wald[k]:7.3f}   likelihood ratio {lr[k]:7.3f}")
```

When the log-likelihood *is* close to quadratic, the three agree closely. For the age
coefficient of @exm-glm-marginals the script obtains the single-coefficient
\( \text{Wald}=2.046 \), \( \text{Score}=2.054 \) and
\( \text{LR}=2.066 \), all far from the \( 5\% \) point \( 3.84 \)
of \( \chi^2(1) \), so the three lead to the same conclusion.

**The score statistic needs only the null fit**, so it can test a term from the model
without it; but it evaluates the information at \( \tilde{\bbeta} \), the wrong place when
the alternative is far away. The recommendation is simple: report likelihood ratio tests,
and where an interval is wanted invert them,
\( \{\beta_j:\text{LR}(\beta_j)\le\chi^2_{1-\alpha}(1)\} \), rather than using the Wald
interval.

## The deviance

::: {#def-glm-deviance}
[Deviance and scaled deviance]

Let \( \hat{\bmu} \) be the fitted means of a model and let \( \theta(\cdot) \) be as
in @def-glm-variance. Assuming every \( y_i\in\mathcal M \), the **saturated model** is the one
that sets \( \mu_i=y_i \). The **deviance** of the fit is
\[
\begin{aligned}
D(\y,\hat{\bmu})&=\sum_{i=1}^n d_i,\\[2pt]
d_i&=2w_i\Bigl[y_i\{\theta(y_i)-\theta(\hat\mu_i)\}-b\{\theta(y_i)\}+b\{\theta(\hat\mu_i)\}\Bigr],
\end{aligned}
\]{#eq-glm-deviance}

and the **scaled deviance** is \( D^{*}=D/\phi=2\{\ell_{\text{sat}}-\ell(\hat{\bmu})\} \).
:::

The deviance does not involve \( \phi \); the scaled deviance does. For the normal family,
\( \theta=\mu \) and \( b(\theta)=\theta^2/2 \) give \( d_i=w_i(y_i-\hat\mu_i)^2 \): the deviance
is the weighted residual sum of squares, the scaled deviance \( \text{SSE}/\sigma^2 \).

| Family | \( d_i/w_i \) |
|---|---|
| normal | \( (y-\hat\mu)^2 \) |
| binomial (\( y \) a proportion) | \( 2\bigl[y\log\dfrac{y}{\hat\mu}+(1-y)\log\dfrac{1-y}{1-\hat\mu}\bigr] \) |
| Poisson | \( 2\bigl[y\log\dfrac{y}{\hat\mu}-(y-\hat\mu)\bigr] \) |
| gamma | \( 2\bigl[-\log\dfrac{y}{\hat\mu}+\dfrac{y-\hat\mu}{\hat\mu}\bigr] \) |
| inverse Gaussian | \( \dfrac{(y-\hat\mu)^2}{y\,\hat\mu^2} \) |

with the convention \( 0\log0=0 \).

::: {#prp-glm-deviance-positive}
[The deviance is nonnegative]

Every \( d_i\ge0 \), with \( d_i=0 \) if and only if \( \hat\mu_i=y_i \). Consequently
\( D\ge0 \), with equality only for the saturated fit, and \( D \) decreases when a model is
enlarged.
:::

::: {.proof}
Write \( \tilde\theta=\theta(y_i) \), \( \hat\theta=\theta(\hat\mu_i) \), so
\( b'(\tilde\theta)=y_i \). Then
\[
\frac{d_i}{2w_i}=y_i(\tilde\theta-\hat\theta)-b(\tilde\theta)+b(\hat\theta)
=b(\hat\theta)-b(\tilde\theta)-b'(\tilde\theta)(\hat\theta-\tilde\theta)\ \ge0
\]
by convexity of \( b \) (@lem-glm-convex): the right-hand side is the Bregman divergence
of \( b \), the amount by which \( b \) lies above its tangent at \( \tilde\theta \). Strict
convexity makes it strict unless \( \hat\theta=\tilde\theta \), that is \( \hat\mu_i=y_i \).
Enlarging the model can only raise the maximized likelihood, and
\( D=2\phi\{\ell_{\text{sat}}-\ell(\hat{\bmu})\} \).
:::

The other standard discrepancy is the **Pearson statistic**,
\[
X^2=\sum_{i=1}^n\frac{w_i(y_i-\hat\mu_i)^2}{V(\hat\mu_i)} ,
\]{#eq-glm-pearson}

which is \( \sum_i(y_i-\hat\mu_i)^2/\widehat{\Var}(Y_i) \) up to the factor \( \phi \). It equals
the deviance for the normal family and approximates it elsewhere, since a second-order
expansion of \( d_i \) about \( \hat\mu_i=y_i \) gives \( w_i(y_i-\hat\mu_i)^2/V(\hat\mu_i) \).

## Analysis of deviance

::: {#thm-glm-deviance}
[Distribution theory for the deviance]

Assume (G1)–(G4).

::: {.enumerate options="label=(\alph*)"}
1. *(Nested models.)* Let \( \C(\X_0)\subseteq\C(\X) \) with \( \rank(\X)-\rank(\X_0)=q \). Under
   the reduced model,
   \[
   \frac{D_0-D}{\phi}\ \xrightarrow{d}\ \chi^2(q),
   \]
   where \( D_0 \) and \( D \) are the deviances of the two fits. When \( \phi \) is unknown it is
   replaced by an estimate, and the ratio \( \{(D_0-D)/q\}/\hat\phi \) is referred to
   \( F(q,n-p) \), by analogy with @thm-glh-f-test.

2. *(Grouped data.)* Suppose the data consist of \( G \) groups, \( G \) and \( p \) fixed, with
   \( y_i \) an average of \( w_i=m_i \) observations and \( m_i\to\infty \) for every \( i \). Then
   under the fitted model \( D/\phi\to\chi^2(G-p) \) and \( X^2/\phi\to\chi^2(G-p) \), so both
   are goodness-of-fit statistics.

3. *(Ungrouped data.)* If the number of groups grows with the sample size, (b) does not
   apply and the limit is generally not \( \chi^2(n-p) \). For binary responses with the
   canonical link the failure is total: \( D \) is a function of \( \hbeta \) alone, so it
   carries no information about fit beyond what \( \hbeta \) already contains.
:::

:::

::: {.proof}
(a) By @def-glm-deviance the saturated log-likelihood cancels in the difference, so
\( (D_0-D)/\phi=2\{\ell(\hbeta)-\ell(\tilde{\bbeta})\} \), which is the likelihood ratio
statistic of @prp-glm-three-tests.

(b) Regard the data as the \( \sum_im_i \) individual observations rather than the \( G \)
averages, and take the saturated model to be the one with a separate parameter per group, so
that its model matrix is the group indicators and it has \( G \) parameters, fixed as
\( m_i\to\infty \). For that model the quantities @eq-glm-leverage evaluated at the individual
observations are \( q_i=1/m_{g(i)} \), with \( g(i) \) the group of observation \( i \), so (G3)
holds; and the information in group \( g \) is \( m_gb''(\theta_g)/\phi\to\infty \), so (G4)
holds. Part (a) applied to the fitted model against the saturated one, with \( q=G-p \), gives
\( D/\phi\to\chi^2(G-p) \). For \( X^2 \), a third-order expansion of \( d_i \) about
\( \hat\mu_i=y_i \) gives \( w_i(y_i-\hat\mu_i)^2/V(\hat\mu_i)+O(w_i\lvert y_i-\hat\mu_i\rvert^3) \);
with \( w_i=m_i \) and \( y_i-\hat\mu_i=O_p(m_i^{-1/2}) \) the remainder summed over the \( G \)
groups is \( O_p(\max_im_i^{-1/2}) \), so \( X^2-D\to0 \) in probability.

(c) For a binary response with the logit link, \( \mathcal M=(0,1) \) and \( y_i\in\{0,1\} \) is
not in \( \mathcal M \); taking limits in @eq-glm-deviance with \( 0\log0=0 \) gives the saturated
log-likelihood \( 0 \) and
\( D=-2\sum_i\{y_i\log\hat\mu_i+(1-y_i)\log(1-\hat\mu_i)\}=-2\ell(\hbeta) \). Writing
\( \hat\eta_i=\x_{(i)}\T\hbeta \) and using \( \log(1-\hat\mu_i)=-\log(1+e^{\hat\eta_i}) \),
\[
\ell(\hbeta)=\sum_i\bigl\{y_i\hat\eta_i-\log(1+e^{\hat\eta_i})\bigr\}
=\hbeta\T\X\T\y-\sum_i\log\bigl(1+e^{\x_{(i)}\T\hbeta}\bigr).
\]
By @cor-glm-marginals, \( \X\T\y=\X\T\hat{\bmu} \), so
\( \ell(\hbeta)=\hbeta\T\X\T\hat{\bmu}-\sum_i\log(1+e^{\x_{(i)}\T\hbeta}) \) depends on the data
only through \( \hbeta \), and two data sets with the same \( \hbeta \) have the same deviance.
:::

::: {#exm-glm-deviance-null}
[Grouped and ungrouped, simulated]

A logistic model with two parameters was fitted to \( 2000 \) data sets
simulated from it, in two formats with the same number of Bernoulli trials over the same
range of \( x \): as \( 40 \) groups of \( 20 \) trials, and as \( 800 \)
single trials. (The listing below repeats the experiment with \( 200 \) replicates, so that
it runs in a few seconds.)

For the grouped data \( G-p=38 \), and the simulated deviances have
mean \( 39.74 \) and standard
deviation \( 9.43 \), against the \( \chi^2(38) \) values \( 38 \) and
\( \sqrt{76}=8.72 \): usable. For the ungrouped data \( n-p=798 \), and the simulated
deviances have mean \( 987.95 \) and standard
deviation \( 20.62 \), against the \( \chi^2(798) \)
values \( 798 \) and \( 39.95 \) — the mean off by nearly
\( 200 \), the spread by a factor of two
([Figure 34.5.2](05-deviance-and-inference.html#fig-glm-deviance-null)). Referring that
deviance to \( \chi^2(n-p) \) would reject a correct model with overwhelming confidence.

Nothing is wrong with the sample size. The two designs carry very nearly the same Fisher
information — the script checks that the standard errors of \( \hbeta \) agree to within
\( 2\% \) — so the estimates behave alike. What fails is the reference distribution, because
the saturated model in the second format has \( 800 \) parameters and grows with the data.
:::

::: {when-format="html"}
![**Figure 34.5.2.** Simulated null distributions of the residual deviance, with the
\( \chi^2 \) density that the naive degrees-of-freedom count suggests. (a) Grouped data: the
approximation works. (b) The same information as ungrouped binary responses: it does
not.](deviance_null.svg){#fig-glm-deviance-null width=100%}
:::

::: {when-format="pdf"}
![Simulated null distributions of the residual deviance, with the
\( \chi^2 \) density that the naive degrees-of-freedom count suggests. (a) Grouped data: the
approximation works. (b) The same information as ungrouped binary responses: it does
not.](deviance_null.pdf){width=100%}
:::

```{.python .run #cell-inference-deviance-null}
def logit_irls(X, y, w, steps=25):
    """Maximum likelihood for a binomial model with logit link; y is a proportion."""
    mu = (w * y + 0.5) / (w + 1.0)
    eta = np.log(mu / (1 - mu))
    beta = np.zeros(X.shape[1])
    for _ in range(steps):
        d = mu * (1 - mu)
        W = w * d
        z = eta + (y - mu) / d
        beta = np.linalg.solve((X.T * W) @ X, (X.T * W) @ z)
        eta = X @ beta
        mu = 1 / (1 + np.exp(-eta))
    return beta, mu


def binomial_deviance(y, mu, w):
    t1 = np.where(y > 0, y * np.log(np.where(y > 0, y, 1.0) / mu), 0.0)
    t2 = np.where(y < 1, (1 - y) * np.log(np.where(y < 1, 1 - y, 1.0) / (1 - mu)), 0.0)
    return 2 * np.sum(w * (t1 + t2))


def null_deviances(X, beta, m, reps, rng):
    """Residual deviance of the fitted model, on `reps` data sets simulated from it."""
    p = 1 / (1 + np.exp(-(X @ beta)))
    w = np.full(len(p), float(m))
    out = np.empty(reps)
    for r in range(reps):
        y = rng.binomial(m, p).astype(float) / m
        _, mu = logit_irls(X, y, w)
        out[r] = binomial_deviance(y, mu, w)
    return out

rng = np.random.default_rng(3407)
beta0 = np.array([0.3, 0.9])
G, m = 40, 20                                       # grouped: 40 groups of 20 trials
Xg = np.column_stack([np.ones(G), np.linspace(-1.5, 1.5, G)])
N = 800                                             # ungrouped: 800 single trials
Xu = np.column_stack([np.ones(N), np.linspace(-1.5, 1.5, N)])

reps = 200
Dg = null_deviances(Xg, beta0, m, reps, rng)
Du = null_deviances(Xu, beta0, 1, reps, rng)
print(f"grouped   : mean {Dg.mean():8.2f} (df {G - 2}), sd {Dg.std():6.2f}"
      f" (sqrt(2 df) {np.sqrt(2 * (G - 2)):.2f})")
print(f"ungrouped : mean {Du.mean():8.2f} (df {N - 2}), sd {Du.std():6.2f}"
      f" (sqrt(2 df) {np.sqrt(2 * (N - 2)):.2f})")
```

::: {.idea}
*Differences* of deviances between nested models are reliable, and are what the
analysis-of-deviance table uses. The deviance *itself* is a goodness-of-fit statistic only
for grouped data with large groups. The two claims are constantly confused because software
prints them in the same table.
:::

The analysis-of-deviance table is built like the sequential analysis-of-variance table
of @def-ss-sequential: add terms in a chosen order, recording the drop in deviance at each
step with its degrees of freedom. The reference distribution is chi-squared when \( \phi \)
is known and \( F \) when it is estimated.

::: {#exm-glm-anodev}
[Analysis of deviance for Engel's data]

Continue the gamma/log model of @exm-glm-engel, with \( u \) the centred logarithm of income
and \( n=235 \). Fitting the three nested models gives

| model | residual df | deviance | change | \( F \) | \( p \) |
|---|---:|---:|---:|---:|---:|
| \( 1 \) | \( 234 \) | \( 39.01 \) | | | |
| \( 1+u \) | \( 233 \) | \( 4.27 \) | \( 34.73 \) | \( 1975.2 \) | \( <10^{-4} \) |
| \( 1+u+u^2 \) | \( 232 \) | \( 4.18 \) | \( 0.0918 \) | \( 5.22 \) | \( 0.023 \) |

Here \( \phi \) is unknown, so each change is divided by \( \hat\phi=0.01758 \)
and referred to \( F(1,232) \), following @thm-glm-deviance(a). One qualification: the log
link is not the gamma's canonical link, so (G1) does not hold and the calibration rests not
on the proof given but on @prp-glm-noncanonical, which extends the limit to this case
without proving it. Income matters overwhelmingly; the quadratic term is borderline. The
residual deviance
\( 4.18 \) is *not* to be compared with \( \chi^2(232) \): the data are
ungrouped, so @thm-glm-deviance(b) does not apply, and the small value reflects the small
dispersion of the gamma fit, not a suspiciously good one.
:::

```{.python .run #cell-deviance-anodev}
import numpy as np
import statsmodels.api as sm
from scipy import stats

engel = sm.datasets.engel.load_pandas().data
y = engel["foodexp"].to_numpy()
u = np.log(engel["income"].to_numpy())
u = u - u.mean()                                   # centred, so the terms are less collinear
n = len(y)

def gamma_deviance(y, mu):
    """D = 2 sum { -log(y/mu) + (y - mu)/mu } for the gamma family."""
    return 2 * np.sum(-np.log(y / mu) + (y - mu) / mu)

models = {"1": np.ones((n, 1)),
          "1 + u": np.column_stack([np.ones(n), u]),
          "1 + u + u^2": np.column_stack([np.ones(n), u, u ** 2])}
fits, dev = {}, {}
for name, Xm in models.items():
    f = sm.GLM(y, Xm, family=sm.families.Gamma(sm.families.links.Log())).fit()
    fits[name], dev[name] = f, gamma_deviance(y, f.fittedvalues)

full = fits["1 + u + u^2"]
mu_full = full.fittedvalues
phi_pearson = np.sum((y - mu_full) ** 2 / mu_full ** 2) / (n - 3)
phi_deviance = dev["1 + u + u^2"] / (n - 3)

names = list(models)
print(f"{'term added':14s} {'df':>4s} {'deviance':>10s} {'change':>9s} {'F':>8s} {'p':>8s}")
prev = None
for name in names:
    df_res = n - models[name].shape[1]
    row = f"{name:14s} {df_res:4d} {dev[name]:10.4f}"
    if prev is not None:
        change = dev[prev] - dev[name]
        F = change / phi_pearson
        row += f" {change:9.4f} {F:8.2f} {stats.f.sf(F, 1, n - 3):8.4f}"
    print(row)
    prev = name
print(f"dispersion: Pearson {phi_pearson:.5f}   deviance {phi_deviance:.5f}")
```

## Estimating the dispersion

For the binomial and Poisson families \( \phi=1 \) by construction; if the data contradict
that, the remedy is a different model, which is Chapter 38. For the normal, gamma and
inverse Gaussian families \( \phi \) is free. Maximum likelihood would take it from the full
likelihood, which for the gamma involves the digamma function and behaves badly under
misspecification. The standard choice is instead the **Pearson estimate**
\[
\hat\phi=\frac{X^2}{n-p}=\frac1{n-p}\sum_{i=1}^n\frac{w_i(y_i-\hat\mu_i)^2}{V(\hat\mu_i)} .
\]{#eq-glm-phi-hat}

The motivation is the moment identity \( \E\{w_i(Y_i-\mu_i)^2/V(\mu_i)\}=\phi \), exact at the
true mean; replacing \( \mu_i \) by \( \hat\mu_i \) costs \( p \) degrees of freedom, by analogy
with \( s^2=\text{SSE}/(n-p) \), the exact normal-theory case (@thm-lm-sigma2). The
alternative \( D/(n-p) \) is also in use; for Engel's data the two give
\( 0.01758 \) and \( 0.01803 \). They can differ
sharply for skewed data, and the Pearson version is the more reliable, since it does not
depend on the far tail of the assumed density.

Because the likelihood equations do not involve \( \phi \) (@thm-glm-score(a)), estimating it
does not disturb \( \hbeta \). It does enter every standard error, through
\( \hat{\A}_n^{-1}=\hat\phi(\X\T\hat{\W}\X)^{-1} \), and its own variability is why the
tests are then referred to an \( F \) distribution.

## Residuals and diagnostics

The raw residual \( y_i-\hat\mu_i \) is uninformative on its own, because its variance
depends on \( \hat\mu_i \). Three standardizations are in use.

::: {#def-glm-residuals}
[Pearson, deviance and Anscombe residuals]

For a fitted generalized linear model, with \( d_i \) the deviance contributions
of @eq-glm-deviance and \( h_{ii} \) the diagonal of the weighted hat matrix
\( \bH=\hat{\W}^{1/2}\X(\X\T\hat{\W}\X)^{-1}\X\T\hat{\W}^{1/2} \), define
\[
\begin{aligned}
r_i^{P}&=\frac{\sqrt{w_i}\,(y_i-\hat\mu_i)}{\sqrt{\hat\phi\,V(\hat\mu_i)}},\qquad
r_i^{D}=\operatorname{sign}(y_i-\hat\mu_i)\sqrt{d_i/\hat\phi},\\[2pt]
r_i^{A}&=\frac{\sqrt{w_i}\,\{A(y_i)-A(\hat\mu_i)\}}{\sqrt{\hat\phi}\;A'(\hat\mu_i)\sqrt{V(\hat\mu_i)}},
\end{aligned}
\]
where \( A(\mu)=\int^{\mu}V(t)^{-1/3}\,dt \). The **standardized** versions divide by
\( \sqrt{1-h_{ii}} \).
:::

The Pearson residuals have squares adding to \( X^2/\hat\phi \), the deviance residuals
squares adding to \( D/\hat\phi \); each splits its own discrepancy measure into
per-observation pieces, so a large deviance residual marks an observation the likelihood
explains badly. The Anscombe residual first applies the transformation that brings the
family closest to normal, \( A(\mu)=\tfrac32\mu^{2/3} \) for the Poisson and
\( 3\mu^{1/3} \) for the gamma.

::: {#exm-glm-residual-shapes}
[Which residual is closest to normal]

Two thousand observations were simulated from a Poisson log-linear model with fitted means
between about \( 1 \) and \( 6 \), and the model refitted.
[Figure 34.5.3](05-deviance-and-inference.html#fig-glm-residuals) shows normal quantile plots
of the three residual types. The Pearson residuals are visibly right-skewed, with sample
skewness \( 0.618 \) — the skewness \( \mu^{-1/2} \)
of @exr-glm-skewness, averaged over the fitted means. The deviance residuals have skewness
\( -0.010 \) and the Anscombe residuals
\( -0.081 \); both are far closer to normal, and neither is uniformly
better. Discreteness remains visible in all three as short horizontal runs. The leverages
sum to \( 2 \) (@prp-proj-leverage(b)) and the largest
is \( 0.0029 \).
:::

::: {when-format="html"}
![**Figure 34.5.3.** Normal quantile plots of the three residual types for a Poisson
log-linear fit with small means. The reference line is the identity. Pearson residuals
inherit the skewness of the family; deviance and Anscombe residuals largely remove
it.](poisson_residuals.svg){#fig-glm-residuals width=100%}
:::

::: {when-format="pdf"}
![Normal quantile plots of the three residual types for a Poisson
log-linear fit with small means. The reference line is the identity. Pearson residuals
inherit the skewness of the family; deviance and Anscombe residuals largely remove
it.](poisson_residuals.pdf){width=100%}
:::

```{.python .run #cell-deviance-residuals}
rng = np.random.default_rng(3405)
N = 2000
Xp = np.column_stack([np.ones(N), rng.uniform(-1, 1, N)])
mu_true = np.exp(Xp @ np.array([1.0, 0.8]))                # means between about 1 and 6
yp = rng.poisson(mu_true).astype(float)
pois = sm.GLM(yp, Xp, family=sm.families.Poisson()).fit()
mu_p = pois.fittedvalues

pearson = (yp - mu_p) / np.sqrt(mu_p)                                   # r = (y - mu)/sqrt(V)
d_i = 2 * (np.where(yp > 0, yp * np.log(np.maximum(yp, 1e-300) / mu_p), 0.0) - (yp - mu_p))
deviance_r = np.sign(yp - mu_p) * np.sqrt(d_i)                          # signed root of d_i
anscombe = 1.5 * (yp ** (2 / 3) - mu_p ** (2 / 3)) / mu_p ** (1 / 6)    # A(mu) = 1.5 mu^{2/3}

W = mu_p                                                   # working weights for Poisson/log
root = np.sqrt(W)[:, None] * Xp
H = root @ np.linalg.inv(root.T @ root) @ root.T           # the weighted hat matrix
h = np.diag(H)

for name, r in [("Pearson", pearson), ("deviance", deviance_r), ("Anscombe", anscombe)]:
    print(f"{name:9s} sd {r.std(ddof=2):6.3f}  skewness {stats.skew(r):7.3f}")
print(f"leverages: sum {h.sum():.4f}, largest {h.max():.4f}")
```

The rest of [Chapter 20](../ch20-residuals-leverage-influence/index.html) transfers under
the substitution \( \X\mapsto\hat{\W}^{1/2}\X \), which is Pregibon's (1981) extension of
the linear-model diagnostics: the leverages of @prp-proj-leverage, a Cook-type
measure \( C_i=(r_i^{P})^2h_{ii}/\{p(1-h_{ii})^2\} \) (@def-res-cooks), and the one-step
deletion approximations of @thm-res-deletion with the working response in place of \( \y \).
The exact deletion identities do not, because dropping a case changes the weights too. For
ungrouped binary data the residuals take only two values at each \( \hat\mu_i \), so only
smoothed plots carry a signal — a point @prp-bin-fit develops.

## Model selection

Because the fit is by maximum likelihood, @def-sel-aic-bic applies directly, with
\( \text{AIC}=-2\ell(\hbeta)+2d \) and \( d=p \) when \( \phi \) is known, \( p+1 \) when it is
estimated; for a \( \phi=1 \) family comparing AIC is comparing \( D+2p \). The rest
of [Chapter 29](../ch29-model-selection/index.html) transfers unchanged,
@prp-sel-selection-bias included.

## Exercises

### A. Check your understanding

::: {#exr-glm-deviance-normal}
[A1]

Verify the normal and Poisson rows of the deviance table from @eq-glm-deviance. For the
Poisson, show that if the model has an intercept then
\( \sum_i(y_i-\hat\mu_i)=0 \) by @cor-glm-marginals, so the deviance reduces to
\( 2\sum_iy_i\log(y_i/\hat\mu_i) \).
:::

::: {#exr-glm-df}
[A2]

A binomial model is fitted to \( 60 \) groups with \( 4 \) parameters, and the residual
deviance is \( 71.2 \) on \( 56 \) degrees of freedom. When is comparing it with
\( \chi^2(56) \) legitimate, and what would you need to know about the data first?
:::

::: {.solution}
By @thm-glm-deviance(b) the comparison needs the number of groups fixed and the group sizes
large, so one needs the \( m_i \): if they are all in the hundreds the comparison is
reasonable, and if many are \( 1 \) or \( 2 \) it is not, whatever the total sample size. A rough
rule is that the expected counts \( m_i\hat\mu_i \) and \( m_i(1-\hat\mu_i) \) should all exceed
about five.
:::

### B. Practice

::: {#exr-glm-wald-invariance}
[B1]

Let \( \hat\beta \) have asymptotic standard error \( s \). Compute the Wald statistic for
\( H_0:\beta=0 \) and, using the delta method, the Wald statistic for the equivalent
hypothesis \( H_0:e^{\beta}=1 \). Show that their ratio is
\( \{(e^{\hat\beta}-1)/\hat\beta\}^2e^{-2\hat\beta} \), and evaluate it at \( \hat\beta=1 \) and
\( \hat\beta=3 \). Why is the likelihood ratio statistic unaffected?
:::

::: {#exr-glm-phi-moment}
[B2]

Show that \( \E\{w_i(Y_i-\mu_i)^2/V(\mu_i)\}=\phi \) exactly, for every observation in an
exponential dispersion family. Deduce that \( X^2/n \) is unbiased for \( \phi \) when the true
means are used, and explain why \( X^2/(n-p) \) is only approximately unbiased when they are
estimated.
:::

::: {#exr-glm-anscombe-poisson}
[B3]

Derive the Anscombe residual for the Poisson family from
\( A(\mu)=\int^{\mu}t^{-1/3}dt \), and check that the standardizing denominator is
\( \hat\mu^{1/6} \). Explain why \( A \) is *not* the variance-stabilizing transformation
\( 2\sqrt\mu \) of @cor-tr-classical, and what each one is optimizing.
:::

::: {.solution}
\( A(\mu)=\tfrac32\mu^{2/3} \), so \( A'(\mu)=\mu^{-1/3} \) and
\( A'(\hat\mu)\sqrt{V(\hat\mu)}=\hat\mu^{-1/3}\hat\mu^{1/2}=\hat\mu^{1/6} \), giving
\( r^A=\tfrac32(y^{2/3}-\hat\mu^{2/3})/\hat\mu^{1/6} \). The variance-stabilizing
transformation solves \( A'\propto V^{-1/2} \) and makes the *variance* constant; the Anscombe
transformation solves \( A'\propto V^{-1/3} \) and makes the *third cumulant*, hence the
skewness, vanish to the order considered. Different objectives, different exponents; the
Anscombe residual is then standardized by hand, so it does not need the variance to be
stable by itself.
:::

::: {#exr-glm-pseudo-r2}
[B4]

Two summaries of predictive power for a generalized linear model are the sample correlation
between \( \y \) and \( \hat{\bmu} \), and the likelihood-based
\( R^2_L=(\ell_M-\ell_0)/(\ell_{\text{sat}}-\ell_0) \), where \( \ell_0 \), \( \ell_M \) and
\( \ell_{\text{sat}} \) are the maximized log-likelihoods of the intercept-only model, the
fitted model and the saturated model. Show that for a family with \( \phi=1 \) the second is
\( (D_0-D)/D_0 \), with \( D_0 \) the null deviance, and deduce that it lies in \( [0,1] \) and
never decreases when a term is added. Why does that last property make it a poor basis for
choosing a model?
:::

::: {.solution}
With \( \phi=1 \), @def-glm-deviance gives \( D=2(\ell_{\text{sat}}-\ell_M) \) and
\( D_0=2(\ell_{\text{sat}}-\ell_0) \), so \( \ell_M-\ell_0=(D_0-D)/2 \) and
\( R^2_L=(D_0-D)/D_0 \). Both deviances are nonnegative by @prp-glm-deviance-positive and
\( D\le D_0 \), since enlarging a model cannot lower the maximized likelihood; so
\( R^2_L\in[0,1] \) and cannot decrease when a term is added. That is the defect of the
unadjusted \( R^2 \) of
[Section 9.5](../ch09-sums-of-squares/05-r-squared.html), the point
@prp-ss-adjusted-r2 makes there: it rewards complexity unconditionally. AIC instead charges
\( 2 \) per parameter against \( -2\ell \), and so can rise.
:::

### C. Going deeper

::: {#exr-glm-deviance-not-fit}
[C1]

Extend @thm-glm-deviance(c) to a Poisson log-linear model in which every \( y_i \) is \( 0 \) or
\( 1 \): show that the deviance is again a function of \( \hbeta \) alone. Then show that the
argument breaks down as soon as some \( y_i\ge2 \), and explain in words what extra information
the larger counts supply. (This is the reason the deviance is a serviceable goodness-of-fit
statistic for Poisson data with large means and a useless one for sparse counts.)
:::

::: {.solution}
For the Poisson deviance,
\( D=2\sum_i\{y_i\log(y_i/\hat\mu_i)-(y_i-\hat\mu_i)\} \). When every \( y_i\in\{0,1\} \) the term
\( y_i\log y_i \) is zero, so \( D=-2\sum_iy_i\log\hat\mu_i+2\sum_i(y_i-\hat\mu_i) \); the second
sum vanishes by @cor-glm-marginals if the model has an intercept, and the first is
\( -2\hbeta\T\X\T\y=-2\hbeta\T\X\T\hat{\bmu} \), again by @cor-glm-marginals. So \( D \) is a
function of \( \hbeta \). If some \( y_i\ge2 \), the term \( \sum_iy_i\log y_i \) is a nonzero
function of the data that is *not* determined by \( \X\T\y \): two data sets with the same
sufficient statistic, one with counts \( (2,0) \) at two identical covariate values and one
with \( (1,1) \), have the same \( \hbeta \) but different deviances. The larger counts carry
information about the spread of the responses within a covariate pattern, which is precisely
what a goodness-of-fit statistic needs and what binary data cannot supply.
:::

