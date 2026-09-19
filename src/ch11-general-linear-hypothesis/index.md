# Testing Models: The General Linear Hypothesis

Almost every question put to a linear model asks whether a simpler model would do: whether a
group of regressors can be dropped, whether two coefficients are equal, whether regions
differ once other variables are accounted for. Each such question names a subspace of the
model space. The whole theory of testing it rests on one picture. The part of the data that
the simpler model cannot explain splits orthogonally into a piece in the *test space* and a
piece in the residual space. The \( F \) statistic compares the squared lengths of the two
pieces, per dimension. Normality turns the comparison into an exact test. This chapter
derives that test in the three forms met in practice, as a comparison of fitted models, as a
statement about estimable functions \( \bLambda\T\bbeta=\bm d \), and as a likelihood ratio.
It shows that the Wald and score tests are the same test in disguise, and proves that no
invariant test is more powerful. It then studies power, which depends on the unknown
parameters only through one noncentrality, and ends with the \( t \) test, the overall
\( F \) test and tests of groups of coefficients, including the ways a group test and the
individual tests within it can disagree. Warnings about what a nonsignificant \( F \) does and
does not show run through the chapter.

**What you need.** The normal linear model and its sampling distributions from
[Chapter 7](../ch07-optimality/index.html) (@eq-opt-normal-model, @thm-opt-sampling, @cor-opt-t, @thm-opt-mle). The noncentral \( \chi^2 \), \( F \) and \( t \) distributions and the \( F \)
statistic for nested subspaces from [Chapter 4](../ch04-quadratic-forms/index.html) (@def-qf-noncentral-f, @thm-qf-ncchisq, @thm-qf-f-power, @thm-qf-nested-f). Nested
projections and the constraint space from [Chapter 6](../ch06-projections/index.html) (@thm-proj-nested, @thm-proj-constraint-space, @thm-proj-M-formula). Estimable functions from
[Chapter 8](../ch08-estimability/index.html) (@def-est-estimable, @thm-est-characterization).
Extra sums of squares and the sum of squares of a linear hypothesis from
[Chapter 9](../ch09-sums-of-squares/index.html) (@def-ss-extra, @thm-ss-hypothesis,
@thm-ss-restricted, @prp-ss-effective-df).

## Roadmap

- [Testing reduced models](01-reduced-models.html)
- [The F statistic as a comparison of projections](02-comparing-projections.html)
- [Testable hypotheses](03-testable-hypotheses.html)
- [The likelihood ratio test and its relatives](04-likelihood-ratio.html)
- [Power and noncentrality](05-power.html)
- [Testing single coefficients and groups of coefficients](06-coefficients.html)
- [Summary and notes](07-summary-and-notes.html)
