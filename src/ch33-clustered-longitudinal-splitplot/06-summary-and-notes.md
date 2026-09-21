# Summary and Notes

## Summary

::: {.idea}

1. A split-plot design randomizes twice and has two error terms: the response
   space splits into a whole-plot stratum with error variance
   \( \sigma_s^2+b\sigma_w^2 \) and a subplot stratum with error variance
   \( \sigma_s^2 \) (@lem-cls-strata).

2. In the balanced design ordinary least squares is still best, because
   \( \C(\V\X)\subseteq\C(\X) \). The whole-plot factor is tested against the
   block-by-treatment mean square on \( (r-1)(a-1) \) degrees of freedom, the
   subplot factor and the interaction against the residual; a comparison of
   whole-plot treatments at one subplot level belongs to neither stratum and needs
   both mean squares (@thm-cls-splitplot). Pooling inflates the whole-plot
   \( F \) and deflates the subplot one.

3. With compound symmetry inside clusters of size \( m \), a cluster-level
   coefficient has its variance multiplied by the design effect
   \( 1+(m-1)\rho \) and a within-cluster coefficient by \( 1-\rho \); the usual
   standard error is wrong by \( \{\E(s^2)/\lambda\}^{1/2} \), with \( \E(s^2) \)
   given exactly (@prp-cls-compound). Clustering moves information from
   between-cluster to within-cluster comparisons.

4. The cluster-robust covariance estimator assumes only that clusters are
   independent, and is consistent as the number of clusters grows, not the number
   of observations. With few clusters it is badly biased downward.

5. The repeated-measures \( F \) test is exact if and only if the covariance
   satisfies sphericity: all pairwise differences of the repeated measurements
   have the same variance (@lem-cls-hf-form, @thm-cls-sphericity). Compound
   symmetry implies sphericity; the converse is false.

6. When sphericity fails, the moment-matched approximation multiplies both
   degrees of freedom by Box's \( \epsilon \) (@prp-cls-box); the
   Greenhouse–Geisser and Huynh–Feldt corrections estimate \( \epsilon \), and
   Hotelling's \( T^2 \) avoids the issue at a cost in
   power (@prp-cls-multivariate).

7. Longitudinal covariance is built from random effects, serial correlation and
   measurement error, and the standard models keep one or two of the
   three (@def-cls-covariance). A wrong covariance model leaves the fixed effects
   unbiased and is repaired by the empirical standard error; restricted likelihoods
   compare only models with the same mean structure (@prp-cls-selection).

8. In the balanced random-coefficient growth model the population curve is the
   average of the subjects' own curves, needing no covariance parameter; the
   predicted subject curve shrinks the subject's own curve towards it by
   \( \bD(\bD+\sigma^2(\X\T\X)^{-1})^{-1} \); and restricted maximum likelihood has
   the closed form of the two-stage moment estimator (@prp-cls-growth).

:::

## Notes and sources

**Split plots.** The design and its two-stratum analysis are due to Yates (1935);
Fisher (1935) had already made the choice of error term the central issue of
experimental design. The route taken here is @lem-cls-strata, which reduces
everything to @thm-tw-balanced plus Kruskal's condition (@thm-proj-kruskal);
Christensen (2020, chapter 11) develops split-plot models from a cluster sampling
model in the same coordinate-free spirit. The cross-stratum simple effect
of @thm-cls-splitplot(f) is standard practical advice given a derivation here.

**Clusters and the design effect.** The design effect \( 1+(m-1)\rho \) is
Kish's (1965), from survey sampling; its damage to standard errors for aggregate
regressors was pressed on econometricians by Moulton (1990). The cluster-robust
estimator @eq-cls-cluster-sandwich is the clustered form of the sandwich
estimators of
[Section 21.4](../ch21-nonnormality-heteroscedasticity-serial/04-sandwich-estimators.html),
whose proof also gives its consistency, cluster by cluster in place of
observation by observation; its small-sample failings and their remedies are
surveyed by Cameron and Miller (2015).

**Sphericity.** Box (1954) computed the effect of an arbitrary covariance on the
two-way analysis of variance and introduced \( \epsilon \) and the approximation
of @prp-cls-box. Huynh and Feldt (1970) isolated the necessary and sufficient
condition of @lem-cls-hf-form, which is named after them, and the fixed-\( N \)
form of @thm-cls-sphericity(c) is theirs; the proof here takes the easier route
through the large-sample limit. Geisser and Greenhouse (1958) showed the
lower-bound reference distribution to be exactly conservative, not merely so
within Box's approximation; Greenhouse and Geisser (1959) made the plug-in and
lower-bound corrections a recommendation, Huynh and Feldt (1976) gave the
bias-corrected @eq-cls-hf-epsilon, and Mauchly (1940) the likelihood ratio test of
sphericity itself. The multivariate alternative is Hotelling's
\( T^2 \) (Hotelling 1931), whose distribution theory we quote from
Anderson (2003).

**Longitudinal covariance.** The three-source decomposition @eq-cls-three-sources
and the empirical variogram are due to Diggle (1988); Diggle, Heagerty, Liang and
Zeger (2002) is the standard account. That a working covariance can be wrong
without invalidating the estimate, provided the standard error is computed
empirically, is from Liang and Zeger (1986), whose generalized estimating
equations are the subject of [Chapter 40](../ch40-glmm-gee/index.html). The rule of @prp-cls-selection(c) follows
from @def-mix-reml and is a common source of error in practice. Degrees of freedom
under an estimated covariance belong to Satterthwaite (1946) and Kenward and
Roger (1997); missing data to [Chapter 41](../ch41-missing-data/index.html), the missing-at-random condition being
Rubin's (1976).

**Growth curves.** Potthoff and Roy (1964) introduced @eq-cls-potthoff-roy and
reduced it to an ordinary multivariate analysis of variance by a weight matrix
fixed in advance; Rao (1965) gave the efficient estimate and recast growth curves
in terms of random coefficients, @eq-cls-random-coefficient; Laird and Ware (1982)
put the two together in the form that became standard software. The shrinkage
@eq-cls-growth-blup is the multivariate signal-to-noise weighting
of [Chapter 27](../ch27-shrinkage/index.html); the closed form
@eq-cls-growth-reml appears to be folklore, proved here by factorizing the
error-contrast likelihood into within- and between-subject pieces.

**The data.** The investment panel of @exm-cls-grunfeld-covariance and
@exm-cls-grunfeld-growth is Grunfeld's (1958), in the eleven-firm version
reconstructed by Kleiber and Zeileis (2010) and distributed with statsmodels — the
panel already fitted with a random firm effect
in [Chapter 32](../ch32-linear-mixed-models/index.html), seen here from two more
angles. The split-plot and repeated-measures examples are simulated.

## References

- Anderson, Theodore W. (2003). *An Introduction to Multivariate Statistical Analysis*. 3rd edition. Hoboken, NJ: Wiley.
- Box, George E. P. (1954). Some Theorems on Quadratic Forms Applied in the Study of Analysis of Variance Problems, II. Effects of Inequality of Variance and of Correlation Between Errors in the Two-Way Classification. *The Annals of Mathematical Statistics* 25(3), 484–498.
- Cameron, A. Colin and Miller, Douglas L. (2015). A Practitioner's Guide to Cluster-Robust Inference. *Journal of Human Resources* 50(2), 317–372.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Diggle, Peter J. (1988). An Approach to the Analysis of Repeated Measurements. *Biometrics* 44(4), 959–971.
- Diggle, Peter J., Heagerty, Patrick, Liang, Kung-Yee and Zeger, Scott L. (2002). *Analysis of Longitudinal Data*. 2nd edition. Oxford: Oxford University Press.
- Fisher, Ronald A. (1935). *The Design of Experiments*. Edinburgh: Oliver and Boyd.
- Geisser, Seymour and Greenhouse, Samuel W. (1958). An Extension of Box's Results on the Use of the \( F \) Distribution in Multivariate Analysis. *The Annals of Mathematical Statistics* 29(3), 885–891.
- Greenhouse, Samuel W. and Geisser, Seymour (1959). On Methods in the Analysis of Profile Data. *Psychometrika* 24(2), 95–112.
- Grunfeld, Yehuda (1958). *The Determinants of Corporate Investment*. PhD thesis, University of Chicago.
- Hotelling, Harold (1931). The Generalization of Student's Ratio. *The Annals of Mathematical Statistics* 2(3), 360–378.
- Huynh, Huynh and Feldt, Leonard S. (1970). Conditions under Which Mean Square Ratios in Repeated Measurements Designs Have Exact \( F \)-Distributions. *Journal of the American Statistical Association* 65(332), 1582–1589.
- Huynh, Huynh and Feldt, Leonard S. (1976). Estimation of the Box Correction for Degrees of Freedom from Sample Data in Randomized Block and Split-Plot Designs. *Journal of Educational Statistics* 1(1), 69–82.
- Kenward, Michael G. and Roger, James H. (1997). Small Sample Inference for Fixed Effects from Restricted Maximum Likelihood. *Biometrics* 53(3), 983–997.
- Kish, Leslie (1965). *Survey Sampling*. New York: Wiley.
- Kleiber, Christian and Zeileis, Achim (2010). The Grunfeld Data at 50. *German Economic Review* 11(4), 404–417.
- Laird, Nan M. and Ware, James H. (1982). Random-Effects Models for Longitudinal Data. *Biometrics* 38(4), 963–974.
- Liang, Kung-Yee and Zeger, Scott L. (1986). Longitudinal Data Analysis Using Generalized Linear Models. *Biometrika* 73(1), 13–22.
- Mauchly, John W. (1940). Significance Test for Sphericity of a Normal \( n \)-Variate Distribution. *The Annals of Mathematical Statistics* 11(2), 204–209.
- Moulton, Brent R. (1990). An Illustration of a Pitfall in Estimating the Effects of Aggregate Variables on Micro Units. *The Review of Economics and Statistics* 72(2), 334–338.
- Potthoff, Richard F. and Roy, Samarendra N. (1964). A Generalized Multivariate Analysis of Variance Model Useful Especially for Growth Curve Problems. *Biometrika* 51(3/4), 313–326.
- Rao, C. Radhakrishna (1965). The Theory of Least Squares When the Parameters Are Stochastic and Its Application to the Analysis of Growth Curves. *Biometrika* 52(3/4), 447–458.
- Rubin, Donald B. (1976). Inference and Missing Data. *Biometrika* 63(3), 581–592.
- Satterthwaite, Franklin E. (1946). An Approximate Distribution of Estimates of Variance Components. *Biometrics Bulletin* 2(6), 110–114.
- Yates, Frank (1935). Complex Experiments. *Supplement to the Journal of the Royal Statistical Society* 2(2), 181–247.
