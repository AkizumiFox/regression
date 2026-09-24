# Multiplicity: Scheffé, Tukey, Bonferroni and the False Discovery Rate

A test with error rate \( 0.05 \) keeps that promise only when it is the only inference drawn from
the data. A regression analysis rarely stops at one: it compares every pair of groups, examines each
of many coefficients, and follows a significant \( F \) test with the contrast that looks most
interesting. This chapter measures the error of a whole family of inferences and controls it in four
ways. Scheffé's method covers every estimable function in a subspace, even functions chosen after
seeing the data, and is exactly the \( F \) test read one direction at a time. Tukey's studentized
range is exact for all pairwise differences of balanced means. Bonferroni's inequality, and Holm's
sharper step-down version, need only the level of each test. The Benjamini–Hochberg procedure
controls the expected proportion of false discoveries, which is what screening many hypotheses
needs. The chapter ends with a guide to choosing among them.

**What you need.** The normal linear model and the sampling distributions of
[Chapter 7](../ch07-optimality/index.html) (@thm-opt-sampling, @cor-opt-t, @cor-opt-quadratic),
independence of projections from
[Chapter 4](../ch04-quadratic-forms/index.html) (@thm-qf-orthogonal-projections), and
orthogonal projections onto subspaces from [Chapter 6](../ch06-projections/index.html) (@thm-proj-nested, @thm-proj-constraint-space).
Estimable functions and contrasts in the
one-way layout come from [Chapter 8](../ch08-estimability/index.html) (@thm-est-oneway-contrasts, @prp-est-contrast-ss),
and the sum of squares of a linear
hypothesis from [Chapter 9](../ch09-sums-of-squares/index.html) (@thm-ss-hypothesis).
The \( F \) test of [Chapter 11](../ch11-general-linear-hypothesis/index.html) (@thm-glh-f-test, @thm-glh-general-f)
and the confidence intervals and ellipsoids of
[Chapter 12](../ch12-intervals-and-bands/index.html) (@thm-ci-estimable-interval, @thm-ci-ellipsoid)
are the single-inference versions of what this chapter does for families. Two pieces of
machinery come in directly: Scheffé's theorem in [Section 13.2](02-scheffe.html) is the
Cauchy–Schwarz inequality together with the generalized Rayleigh quotient of
[Chapter 1](../ch01-matrix-algebra/index.html) (@prp-mat-cauchy-schwarz,
@cor-mat-generalized-rayleigh), and [Section 13.4](04-bonferroni-holm.html) uses the
independence of uncorrelated normal blocks from
[Chapter 3](../ch03-multivariate-normal/index.html) (@thm-mvn-independence). The independence of
orthogonal projections (@thm-qf-orthogonal-projections) is named above for the joint distribution
it justifies, but no proof in the chapter cites it.

## Roadmap

- [The multiple testing problem](01-multiple-testing.html)
- [Scheffé's projection method](02-scheffe.html)
- [Tukey's studentized range](03-tukey.html)
- [Bonferroni and Holm](04-bonferroni-holm.html)
- [False discovery rate control](05-fdr.html)
- [Choosing a procedure](06-choosing.html)
- [Summary and notes](07-summary-and-notes.html)
