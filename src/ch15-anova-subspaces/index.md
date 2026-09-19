# Analysis of Variance as Subspace Comparison

An experimenter assigns treatments to units at random and measures a response. The simplest such
experiment has one factor, and its analysis is the prototype for Part IV. The mean vector lies in the
space of vectors constant within groups, that space splits orthogonally into the constants and a
treatment space, and every question about the treatments concerns the component of the mean in the
treatment space. This chapter turns that geometry into a working analysis. It states the one-way model and what its parameters mean, shows that single contrasts, families of
contrasts and the overall \( F \) test are projections of the data onto pieces of one subspace, and
derives the analysis of variance table and its expected mean squares. It splits the treatment sum of squares into single degrees of freedom along directions
chosen for their meaning, including polynomial trends in a quantitative factor. Finally it asks what
balance buys: plain orthogonality, the largest guaranteed power, and a test that tolerates unequal
variances. When the design is unbalanced and the variances differ, the classical test can go badly
wrong, and the chapter describes the heteroscedasticity-robust tests of Welch and of Brown and
Forsythe.

**What you need.** The one-way projection and its generalized inverses from
[Chapter 6](../ch06-projections/index.html) (@exm-proj-oneway-rank, @exm-proj-oneway-M, @thm-proj-nested);
estimability, side conditions and contrasts from
[Chapter 8](../ch08-estimability/index.html) (@thm-est-oneway-contrasts, @prp-est-contrast-ss, @thm-est-side-conditions);
orthogonal decompositions and expected mean squares from
[Chapter 9](../ch09-sums-of-squares/index.html) (@thm-ss-expected-mean-squares, @thm-ss-single-df);
the \( F \) test and its power from [Chapter 11](../ch11-general-linear-hypothesis/index.html) (@thm-glh-f-test, @thm-glh-power);
the \( t \) interval from
[Chapter 12](../ch12-intervals-and-bands/index.html) (@thm-ci-estimable-interval); and the
simultaneous methods of [Chapter 13](../ch13-multiplicity/index.html) (@thm-mc-scheffe, @thm-mc-tukey, @thm-mc-bonferroni).
The pure-error lack-of-fit test of
[Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html) (@thm-cor-lack-of-fit) reappears
as trend analysis.

## Roadmap

- [The one-way layout as a linear model](01-oneway-model.html)
- [Estimable contrasts](02-contrasts.html)
- [The analysis of variance table from projections](03-anova-table.html)
- [Orthogonal contrasts](04-orthogonal-contrasts.html)
- [Balanced and unbalanced one-way layouts](05-balance.html)
- [Summary and notes](06-summary-and-notes.html)
