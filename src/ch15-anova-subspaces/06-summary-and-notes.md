# Summary and Notes

## Summary

::: {.idea}

1. The one-way layout has mean space \( \C(\Z) \), the vectors constant within groups, in both its forms (@def-aov-oneway). The fit is
   the group means, and \( s^2 \) pools the group variances. Different side conditions define different "overall means" when the layout
   is unbalanced ([Section 15.1](01-oneway-model.html)).

2. Contrasts are the directions of the treatment space \( \mathcal T=\C(\Z)\cap\bone\perpc \), and a contrast's test, interval and sum of
   squares come from the projection onto its line (@prp-aov-contrast-ss). Families use Bonferroni or Dunnett for comparisons with a
   control, Tukey for all pairs and Scheffé for contrasts suggested by the data ([Section 15.2](02-contrasts.html)).

3. A contrast test captures the fraction \( \cos^2\varphi \) of the overall noncentrality, \( \varphi \) being its angle with the true treatment
   component (@prp-aov-contrast-angle): a well-chosen planned contrast beats the overall test, a badly chosen one has no power.

4. The decomposition \( \I=\bP_0+(\M-\bP_0)+(\I-\M) \) gives the analysis of variance table. The \( F \) statistic has the
   \( F(g-1,n-g,\gamma) \) law with \( \gamma=\sum_kn_k(\mu_k-\bar{\mu})^2/\sigma^2 \), and it is the likelihood ratio test of equal means
   (@thm-aov-oneway-f). Without normality, \( \E(\text{MSE})=\sigma^2 \) and \( \E(\text{MSB})=\sigma^2+\sum_kn_k(\mu_k-\bar{\mu})^2/(g-1) \)
   (@thm-aov-oneway-ems). The quantity \( (g-1)F \) is the largest squared contrast \( t \) statistic.

5. A complete set of orthogonal contrasts, orthogonal in the sense \( \sum_kc_kd_k/n_k=0 \), splits \( \text{SSB} \) into \( g-1 \) independent
   single-degree-of-freedom pieces whose noncentralities add up to \( \gamma \), and \( F \) is the average of their \( F \) statistics
   (@thm-aov-orthogonal-contrasts). Successive splits of the treatments give such sets for any group sizes (@prp-aov-tree-contrasts), and
   polynomial contrasts give a trend analysis, whose higher-degree pieces are the pure-error lack-of-fit test of a low-degree curve
   (@exm-aov-nitrogen-trend).

6. Balance makes design orthogonality plain orthogonality and maximizes the noncentrality guaranteed against every configuration of
   means with a given range (@prp-aov-balance-power).

7. With unequal variances, \( \text{MSB} \) and \( \text{MSE} \) estimate the same average variance under equal means if the layout is balanced,
   but not otherwise. The classical test is liberal when the small groups are the variable ones and conservative when the large groups are
   (@thm-aov-heteroscedastic). The tests of Welch and of Brown and Forsythe hold their level approximately, and Welch–Satterthwaite
   intervals do the same for contrasts ([Section 15.5](05-balance.html)).

8. Residuals have variance \( \sigma^2(1-1/n_k) \); plotted against the group means and on a normal quantile plot after studentizing, they
   check the variances and normality. Within-group correlation, which the data cannot reveal, is the most damaging failure.

:::

## Notes and sources

**Origins.** The analysis of variance is Fisher's. It appears in *Statistical Methods for Research Workers* (Fisher 1925), developed
for the field experiments at Rothamsted, and *The Design of Experiments* (Fisher 1935) made randomization the logical basis of the
analysis. The derivation of the one-way model from unit–treatment additivity and random assignment, sketched in
[Section 15.1](01-oneway-model.html), is developed in detail by Kempthorne (1952). Scheffé (1959) remains the classic treatment of the fixed-effects theory,
including the effects of departures from the assumptions.

**The textbook treatments.** Christensen (2020, chapter 4) derives the one-way analysis entirely from the projections \( \M \) and
\( \M-\bP_0 \), including the characterization of contrasts through the space onto which \( \M-\bP_0 \) projects; this chapter follows his
geometric approach, with the treatment space as the organizing object. Rencher and Schaalje (2008, chapter 13) give the balanced case
through full-and-reduced models and through the general linear hypothesis, derive the expected mean squares both ways, and treat
orthogonal and orthogonal polynomial contrasts, including the result that the overall \( F \) statistic averages the single-degree-of-freedom
statistics of a complete orthogonal set. Their section 12.7 covers the testable hypotheses of non-full-rank models used in @thm-aov-oneway-f.
Seber and Lee (2003, section 8.2) discuss simultaneous intervals and the assumptions, including unequal variances. Note that Rencher and Schaalje and Christensen define
noncentrality as half of our \( \gamma \).

**Contrasts and multiple comparisons.** Tables of orthogonal polynomials for equally spaced levels go back at least to Fisher and
Yates (1938). The exact method for comparisons with a control is due to Dunnett (1955), and the conservativeness of the Tukey–Kramer
intervals for unequal group sizes was proved by Hayter (1984); [Chapter 13](../ch13-multiplicity/index.html) (@exr-mc-control) characterizes the many-to-one multiplier and gives further references. @prp-aov-contrast-angle and @prp-aov-tree-contrasts are elementary consequences of the projection picture. Part (a) of @prp-aov-balance-power is @exr-glh-unbalanced-lf of Chapter 11; part (b), the optimality of balance, is a short calculation.

**Unequal variances.** Box (1954) computed the effect of unequal variances on the one-way \( F \) test, showing that it is small when the
groups are of equal size and can be large otherwise, in the directions described by @thm-aov-heteroscedastic. The test of Welch (1951)
weights the groups by their estimated precisions; the test of Brown and Forsythe (1974a) keeps the unweighted numerator and estimates its
null expectation. Welch obtained his correction factor and degrees of freedom from a series expansion of the statistic's distribution
to first order in \( 1/(n_k-1) \); the Brown–Forsythe degrees of freedom, like those of the contrast intervals, are Satterthwaite's
(1946) moment match. The robust test of equal variances by
analysis of variance of absolute deviations is due to Levene (1960), and the version based on deviations from medians to Brown and
Forsythe (1974b). The simulation of @exm-aov-size-simulation is ours. [Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) returns to heteroscedasticity with weighted least squares and
heteroscedasticity-consistent standard errors, and Chapter 32 to random effects.

## References

- Box, George E. P. (1954). Some Theorems on Quadratic Forms Applied in the Study of Analysis of Variance Problems, I. Effect of Inequality of Variance in the One-Way Classification. *The Annals of Mathematical Statistics* 25(2), 290–302.
- Brown, Morton B. and Forsythe, Alan B. (1974a). The Small Sample Behavior of Some Statistics Which Test the Equality of Several Means. *Technometrics* 16(1), 129–132.
- Brown, Morton B. and Forsythe, Alan B. (1974b). Robust Tests for the Equality of Variances. *Journal of the American Statistical Association* 69(346), 364–367.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Dunnett, Charles W. (1955). A Multiple Comparison Procedure for Comparing Several Treatments with a Control. *Journal of the American Statistical Association* 50(272), 1096–1121.
- Fisher, Ronald A. (1925). *Statistical Methods for Research Workers*. Edinburgh: Oliver and Boyd.
- Fisher, Ronald A. (1935). *The Design of Experiments*. Edinburgh: Oliver and Boyd.
- Fisher, Ronald A. and Yates, Frank (1938). *Statistical Tables for Biological, Agricultural and Medical Research*. Edinburgh: Oliver and Boyd.
- Hayter, Anthony J. (1984). A Proof of the Conjecture That the Tukey–Kramer Multiple Comparisons Procedure Is Conservative. *The Annals of Statistics* 12(1), 61–75.
- Kempthorne, Oscar (1952). *The Design and Analysis of Experiments*. New York: Wiley.
- Levene, Howard (1960). Robust Tests for Equality of Variances. In I. Olkin et al. (eds), *Contributions to Probability and Statistics: Essays in Honor of Harold Hotelling*, 278–292. Stanford, CA: Stanford University Press.
- Rencher, Alvin C. and Schaalje, G. Bruce (2008). *Linear Models in Statistics*. 2nd edition. Hoboken, NJ: Wiley.
- Satterthwaite, Franklin E. (1946). An Approximate Distribution of Estimates of Variance Components. *Biometrics Bulletin* 2(6), 110–114.
- Scheffé, Henry (1959). *The Analysis of Variance*. New York: Wiley.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Welch, Bernard L. (1951). On the Comparison of Several Mean Values: An Alternative Approach. *Biometrika* 38(3/4), 330–336.
