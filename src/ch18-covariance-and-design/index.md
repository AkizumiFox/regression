# Analysis of Covariance and Designed Experiments

An experiment is designed before any response is measured. The experimenter decides which units
receive which treatment, which nuisance variation to remove by grouping units into blocks, and
which extra measurements to take on each unit before treatment begins. Each decision shows up in
the linear model: blocks as additional factors, extra measurements as covariates, and the
randomization as the reason the error assumptions are believable. This chapter works through the
standard designed structures with the projection tools of Parts II and III. Analysis of covariance
turns out to be the Frisch–Waugh–Lovell theorem applied to a factor model. Adjusted means are
estimable functions with a correction term. The test of parallel slopes compares nested subspaces.
Randomized blocks, Latin squares and balanced incomplete blocks are layouts whose balance makes the
projections simple. A lost observation is handled exactly by least squares on the observations that
remain, and the classical ways of filling in the gap are shortcuts to the same answer.

**What you need.** [Chapter 6](../ch06-projections/index.html) (nested projections, @thm-proj-nested, and the Frisch–Waugh–Lovell theorem, @thm-proj-fwl),
[Chapter 8](../ch08-estimability/index.html) (estimability, and the parallel and separate lines models of
[Section 8.6](../ch08-estimability/06-several-factors.html)),
[Chapter 9](../ch09-sums-of-squares/index.html) (extra sums of squares and orthogonal designs),
[Chapter 11](../ch11-general-linear-hypothesis/index.html) and [Chapter 13](../ch13-multiplicity/index.html)
(\( F \) tests and simultaneous inference), [Chapter 12](../ch12-intervals-and-bands/index.html) (\( t \) intervals,
prediction intervals and the Working–Hotelling band), and, from this part,
[Chapter 15](../ch15-anova-subspaces/index.html) and [Chapter 16](../ch16-multiway-layouts/index.html)
(the one-way and additive two-way layouts) and [Chapter 17](../ch17-unbalanced-data/index.html) (cell means and least
squares means in unbalanced layouts). [Section 18.1](01-ancova.html) also uses the distribution of
the sample variance from [Chapter 4](../ch04-quadratic-forms/index.html) (@cor-qf-sample-variance).
Of the chapters named, [Chapter 15](../ch15-anova-subspaces/index.html) is the one the designs here
are read against rather than argued from: nothing in this chapter cites it.

## Roadmap

- [Analysis of covariance as a projection problem](01-ancova.html)
- [Adjusted treatment means](02-adjusted-means.html)
- [Heterogeneous slopes](03-heterogeneous-slopes.html)
- [Randomized complete blocks](04-randomized-blocks.html)
- [Latin squares and incomplete blocks](05-latin-squares-incomplete-blocks.html)
- [Missing observations in designed experiments](06-missing-observations.html)
- [Summary and notes](07-summary-and-notes.html)
