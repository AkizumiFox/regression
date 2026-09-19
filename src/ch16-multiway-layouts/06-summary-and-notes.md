# Summary and Notes

## Summary

::: {.idea}

1. In a balanced layout every projection of interest is a Kronecker product of averaging matrices \( \bar{\mathbf{J}} \), centring matrices \( \mathbf{C} \) and identities. Products that differ by \( \bar{\mathbf{J}} \) against \( \mathbf{C} \) in some position are orthogonal, and ranks multiply (@lem-tw-kron).

2. The additive two-way model has projection \( \bP_0+\bP_A+\bP_B \). Row contrasts are estimated and judged as in a one-way layout that ignores the columns. With replication, the residual splits into interaction (lack of fit) and pure error. With one observation per cell, Tukey's one-degree-of-freedom test is an exact test against the kind of nonadditivity a wrong scale produces, and its slope suggests a power transformation (@thm-tw-additive).

3. Interaction is a property of the table of cell means. It is present iff some interaction contrast (a zero-margin combination of cell means) is nonzero, and it has \( (a-1)(b-1) \) dimensions whatever parameterization is used (@def-tw-interaction, @prp-tw-interaction-equiv).

4. An increasing transformation can remove interaction only if the orders of the levels are the same everywhere. Qualitative (crossing) interaction cannot be removed. Main effects in the presence of interaction are weighted averages of simple effects, and their values, even their signs, depend on the weights (@prp-tw-removable, @prp-tw-weighted-main).

5. The balanced two-way layout with replication splits into grand mean, \( A \), \( B \), \( AB \) and error, with ranks \( 1 \), \( a-1 \), \( b-1 \), \( (a-1)(b-1) \), \( ab(m-1) \). The \( F \) tests are exact, and the sequential, Type II and Type III sums of squares all coincide (@thm-tw-balanced).

6. Products of orthogonal contrasts split the interaction sum of squares into single degrees of freedom. Tukey's method covers the levels of one factor, Scheffé's all interaction contrasts, and Bonferroni or Holm planned families (@prp-tw-product-contrasts).

7. A balanced factorial with \( t \) factors has one term for each subset of factors, with projection \( \bP_S \) built from centring matrices in the positions of \( S \). Hierarchical models have order-free sums of squares, and \( \bP_S\y \) is an alternating sum of marginal means (@thm-tw-factorial). The three-factor interaction measures how a two-factor interaction changes across the third factor (@prp-tw-three-factor).

8. In a \( 2^t \) factorial every term has one degree of freedom, and its effect is a difference of two means with variance \( 4\sigma^2/n \) (@prp-tw-two-level). Pooling terms into error is exact if they are null and chosen in advance, and conservative if they are chosen in advance but not null. It is invalid if they are chosen because they look small (@prp-tw-pooling).

9. A nested factor has no main effect. Its sum of squares is the crossed \( B \) plus \( AB \), the variation among inner levels within each outer level (@def-tw-nested, @thm-tw-nested). When the inner levels are a sample, the outer factor must be tested against the inner mean square, not against the within-level error.

:::

## Notes and sources

**Sources.** The coverage of this chapter was checked against Rencher and Schaalje (2008, chapter 14), Christensen (2020, chapter 7) and Seber and Lee (2003, sections 8.4–8.6). The projection approach follows the vector-space treatment of Christensen (2020), which develops interaction contrasts and their product construction (@prp-tw-product-contrasts) with particular care. Rencher and Schaalje compute the expected mean squares directly, and projections replace that computation here. Scheffé (1959) remains the classical account of fixed-effects multiway layouts, including main effects and interactions defined with general weights (compare @prp-tw-interaction-equiv(v)).

**Factorial experiments.** Fisher (1935) argued for factorial experiments over one-factor-at-a-time trials: every observation contributes to every main effect, and interactions can be studied at all. Yates (1937) developed the analysis of \( 2^t \) and \( 3^t \) designs, including the algorithm for \( 2^t \) effects that carries his name. Box, Hunter and Hunter (2005) and Wu and Hamada (2009) are standard modern texts.

**Tukey's test.** The one-degree-of-freedom test is due to Tukey (1949). Its exact distribution follows from the theorem of Milliken and Graybill (1970) on regressors computed from the fitted values (@thm-cor-fitted-regressors), because Tukey's regressor is, after the additive model is projected out, the table of squared fitted values; @exr-cor-tukey-nonadditivity made the same connection. The power rule @eq-tw-tukey-power is a heuristic, and Chapter 22 treats transformations systematically.

**Interaction and scale.** The dependence of interaction on the scale of measurement, and the order condition of @prp-tw-removable, belong to the theory of conjoint measurement of Luce and Tukey (1964). The counterexample of @exr-tw-nonremovable is an instance of their double cancellation condition. In clinical trials a treatment effect that changes sign across patient subgroups is called a *qualitative interaction*, and Gail and Simon (1985) gave a likelihood ratio test for it.

**Unreplicated factorials and pooling.** The half-normal plot is due to Daniel (1959), and the pseudo standard error of @exr-tw-lenth to Lenth (1989). The effect of preliminary tests on subsequent inference was first studied by Bancroft (1944), and the simulation of [Figure 16.4.2](04-higher-way.html#fig-tw-pooling) illustrates his conclusions in a simple case.

**Nested designs.** Hurlbert (1984) named *pseudoreplication*, treating observations within an experimental unit as independent replicates of the unit's treatment, and documented how widespread it was in ecology. Random inner levels are taken up with variance components in Chapters 32 and 33. Unbalanced multiway data are the subject of [Chapter 17](../ch17-unbalanced-data/index.html), and randomized blocks and Latin squares, additive layouts whose structure comes from randomization, of [Chapter 18](../ch18-covariance-and-design/index.html).

## References

- Bancroft, T. A. (1944). On Biases in Estimation Due to the Use of Preliminary Tests of Significance. *The Annals of Mathematical Statistics* 15(2), 190–204.
- Box, George E. P., Hunter, J. Stuart and Hunter, William G. (2005). *Statistics for Experimenters: Design, Innovation, and Discovery*. 2nd edition. Hoboken, NJ: Wiley.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Daniel, Cuthbert (1959). Use of Half-Normal Plots in Interpreting Factorial Two-Level Experiments. *Technometrics* 1(4), 311–341.
- Fisher, Ronald A. (1935). *The Design of Experiments*. Edinburgh: Oliver and Boyd.
- Gail, Mitchell and Simon, Richard (1985). Testing for Qualitative Interactions between Treatment Effects and Patient Subsets. *Biometrics* 41(2), 361–372.
- Hurlbert, Stuart H. (1984). Pseudoreplication and the Design of Ecological Field Experiments. *Ecological Monographs* 54(2), 187–211.
- Lenth, Russell V. (1989). Quick and Easy Analysis of Unreplicated Factorials. *Technometrics* 31(4), 469–473.
- Luce, R. Duncan and Tukey, John W. (1964). Simultaneous Conjoint Measurement: A New Type of Fundamental Measurement. *Journal of Mathematical Psychology* 1(1), 1–27.
- Milliken, George A. and Graybill, Franklin A. (1970). Extensions of the General Linear Hypothesis Model. *Journal of the American Statistical Association* 65(330), 797–807.
- Rencher, Alvin C. and Schaalje, G. Bruce (2008). *Linear Models in Statistics*. 2nd edition. Hoboken, NJ: Wiley.
- Scheffé, Henry (1959). *The Analysis of Variance*. New York: Wiley.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Tukey, John W. (1949). One Degree of Freedom for Non-Additivity. *Biometrics* 5(3), 232–242.
- Wu, C. F. Jeff and Hamada, Michael S. (2009). *Experiments: Planning, Analysis, and Optimization*. 2nd edition. Hoboken, NJ: Wiley.
- Yates, Frank (1937). *The Design and Analysis of Factorial Experiments*. Technical Communication No. 35. Harpenden: Imperial Bureau of Soil Science.
