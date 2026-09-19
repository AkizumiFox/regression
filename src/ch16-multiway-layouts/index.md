# Multi-Way Layouts and Interaction

An experiment with several factors asks more than one question. It asks whether each factor matters, and also whether the effect of one factor depends on the level of another. In a balanced layout, where every combination of levels is observed equally often, these questions correspond to mutually orthogonal subspaces of observation space. Each subspace has a projection that is a Kronecker product of averaging and centring matrices, and the analysis of variance is Pythagoras applied to that decomposition. This chapter builds the decomposition for two factors, with and without replication, for any number of crossed factors, and for factors nested inside other factors. The algebra is the easy part. The harder part is interpretation: what an interaction is, when a transformation removes it, what a main effect means when interaction is present, and when it is safe to pool one term into another.

**What you need.** Kronecker products ([Section 1.10](../ch01-matrix-algebra/10-kronecker.html)); independence of projections of a normal vector and the noncentral \( F \) distribution ([Chapter 4](../ch04-quadratic-forms/index.html)); nested projections and the two-factor projection identity of @exm-proj-two-factor ([Section 6.5](../ch06-projections/05-nested.html)); estimability with two factors ([Section 8.6](../ch08-estimability/06-several-factors.html)); analysis of variance decompositions, orthogonal designs and expected mean squares ([Chapter 9](../ch09-sums-of-squares/index.html)); \( F \) tests of reduced models ([Chapter 11](../ch11-general-linear-hypothesis/index.html)); simultaneous inference ([Chapter 13](../ch13-multiplicity/index.html)); the lack-of-fit test and regressors built from fitted values ([Sections 14.6](../ch14-correlation-lack-of-fit-prediction/06-lack-of-fit.html) and [14.7](../ch14-correlation-lack-of-fit-prediction/07-near-replicates.html)); and the one-way analysis of variance ([Chapter 15](../ch15-anova-subspaces/index.html)).

## Roadmap

- [Two-way additive models](01-additive.html)
- [Interaction and its interpretation](02-interaction.html)
- [Balanced two-way layouts via orthogonal projections](03-balanced.html)
- [Higher-way layouts](04-higher-way.html)
- [Nested factors](05-nested.html)
- [Summary and notes](06-summary-and-notes.html)
