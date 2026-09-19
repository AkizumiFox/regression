# Summary and notes

## Summary

::: {.idea}

1. A linear hypothesis says that the mean vector lies in a subspace \( \C(\X_0) \) of the model space.
   Its \( F \) statistic compares the extra sum of squares \( \norm{(\M-\Mo)\y}^2 \), per dimension of the
   test space, with \( s^2 \). Under normal errors it is \( F(r-r_0,n-r,\gamma) \), central under the hypothesis,
   with \( \gamma \) the squared distance from the mean to the reduced model in units of \( \sigma^2 \)
   (@def-glh-reduced-model, @thm-glh-f-test, @prp-glh-p-value).

2. Geometrically, \( F \) is a function of the angle between the reduced model's residual vector and the test
   space. The rejection region is a cone, the statistic is invariant under shifts within the reduced model,
   rescaling and rotations of the test and residual spaces, and in canonical coordinates every hypothesis
   is a test that a \( q \)-dimensional normal mean is zero. The null distribution holds for every
   spherically symmetric error law (@prp-glh-angle, @prp-glh-invariance, @prp-glh-canonical,
   @prp-glh-spherical).

3. A constraint \( \bLambda\T\bbeta=\bm d \) is a hypothesis about the data only if it is estimable and consistent.
   Nonestimable constraints cannot be tested at all. For a testable hypothesis,
   \( (\bLambda\T\hbeta-\bm d)\T[\bLambda\T(\X\T\X)\ginv\bLambda]\ginv(\bLambda\T\hbeta-\bm d) \) is the extra sum of squares of
   the reduced model with offset \( \X\bb_0 \), whatever the rank of \( \X \) or \( \bLambda \), and restricted least
   squares has a closed form through any generalized inverse (@def-glh-testable, @prp-glh-nonestimable,
   @thm-glh-general-f, @prp-glh-restricted-ls).

4. The \( F \) test is the likelihood ratio test. The Wald, likelihood ratio and score statistics are
   \( nu \), \( n\log(1+u) \) and \( nu/(1+u) \) with \( u=qF/(n-r) \), so they always satisfy
   Wald \( > \) LR \( > \) Score. Calibrated by \( \chi^2 \), they disagree and have inflated sizes in small samples
   (@thm-glh-lrt, @prp-glh-trinity).

5. The noncentral \( F \) family has monotone likelihood ratio, and every test invariant under the symmetries
   of the problem is a function of \( F \). Hence the \( F \) test is uniformly most powerful among invariant
   tests (@lem-glh-mlr, @thm-glh-ump-invariant).

6. Power depends only on \( \gamma \), increases with it, and grows with replication. In the one-way layout
   the least favourable configuration for a given range puts two means at the extremes and the rest at the
   midpoint. Observed power is a function of the p-value and adds nothing to it
   (@thm-glh-power, @prp-glh-least-favourable, @prp-glh-observed-power).

7. For one estimable function the \( F \) test is the two-sided \( t \) test, and the one-sided version has
   exact size for one-sided hypotheses. A group statistic \( qF \) is the largest squared \( t \) statistic over all
   combinations in the group, which explains why group and individual tests can disagree in both directions
   (@thm-glh-t-test, @prp-glh-max-t).

:::

## Notes and sources

**Origins.** The ratio of mean squares as a test statistic is Fisher's. It grew out of his analysis of
variance in the early 1920s, and Fisher (1925) tabulated it on the logarithmic scale, as
\( z=\tfrac12\log F \). Snedecor (1934) introduced the letter \( F \) in Fisher's honour and published tables of the
ratio itself. The exact \( t \) distribution of a regression coefficient, with \( n-p \) degrees of freedom, is also due to
Fisher (1922). The geometry of [Section 11.2](02-comparing-projections.html) follows the coordinate-free tradition
described in the notes to [Chapter 6](../ch06-projections/index.html). Christensen (2020, chapter 3)
develops testing entirely in that language, with test spaces, offsets and nonestimable constraints. Rencher and Schaalje (2008, chapters 8
and 12) and Seber and Lee (2003, chapter 4) give the matrix treatment for full-rank and rank-deficient models,
including the canonical elimination of @exr-glh-eliminate. Sen and Srivastava (1990, chapter 3) and Agresti
(2015, sections 3.1–3.2) emphasize applied examples. Note that Rencher and Schaalje and Christensen define the noncentrality as half of our \( \gamma \)
([Section 4.1](../ch04-quadratic-forms/01-noncentral.html)).

**Testable hypotheses.** The insistence that only estimable functions can be tested goes back at least to
Searle (1971), whose treatment of "testable hypotheses" in models not of full rank is the standard reference;
@prp-glh-nonestimable is a minimal version of his argument. Christensen (2020, section 3.3) shows how a
nonestimable constraint reduces to its estimable part, which in our treatment is @prp-ss-effective-df.

**Likelihood ratio, Wald and score.** Neyman and Pearson (1933) introduced the most powerful test and the
likelihood ratio principle for composite hypotheses. The large-sample \( \chi^2 \) calibration of the likelihood
ratio is due to Wilks (1938), the Wald statistic to Wald (1943), and the score statistic to Rao (1948). The
inequality Wald \( \ge \) LR \( \ge \) score for linear restrictions in normal regression, and the resulting
possibility of conflicting verdicts, were pointed out by Berndt and Savin (1977) and discussed further by
Breusch (1979). The "\( nR^2 \)" form of the score test in @exr-glh-nr2 is widely used in econometrics.

**Invariance.** Lehmann and Romano (2005) develop the theory of invariant tests and prove that the \( F \)
test is uniformly most powerful invariant; @thm-glh-ump-invariant follows their approach in outline.
Scheffé (1959) treats the power of \( F \) tests in the analysis of variance. The robustness of the null distribution to spherically symmetric errors (@prp-glh-spherical) is part
of the theory of spherical and elliptical distributions surveyed by Fang, Kotz and Ng (1990).

**Power.** The noncentral \( F \) distribution was tabulated for power calculations by Tang (1938). The
least favourable configuration of @prp-glh-least-favourable is the basis of the classical sample-size tables
for the analysis of variance. Cohen (1988) popularized effect-size based power analysis. Hoenig and Heisey (2001) explain why "observed power" computed from the data is uninformative, the
point of @prp-glh-observed-power. Wasserstein and Lazar (2016) discuss common misreadings of p-values. Equivalence tests, which
can support a claim of "no important difference", go back to Schuirmann (1987).

**Other.** The test of equality of two regressions in @exr-glh-chow and @exr-glh-chow-few is due to Chow (1960).
Longley (1967) assembled his data to test the numerical accuracy of regression programs
([Chapter 10](../ch10-computation/index.html)).
Monotonicity properties of the \( \chi^2 \), \( F \) and \( t \) distributions in their degrees of freedom, of the kind behind the dilution effect of @exr-glh-critical-growth, are studied systematically by Ghosh (1973).

## References

- Agresti, Alan (2015). *Foundations of Linear and Generalized Linear Models*. Hoboken, NJ: Wiley.
- Berndt, Ernst R. and Savin, N. Eugene (1977). Conflict among Criteria for Testing Hypotheses in the Multivariate Linear Regression Model. *Econometrica* 45(5), 1263–1277.
- Breusch, Trevor S. (1979). Conflict among Criteria for Testing Hypotheses: Extensions and Comments. *Econometrica* 47(1), 203–207.
- Chow, Gregory C. (1960). Tests of Equality between Sets of Coefficients in Two Linear Regressions. *Econometrica* 28(3), 591–605.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Cohen, Jacob (1988). *Statistical Power Analysis for the Behavioral Sciences*. 2nd edition. Hillsdale, NJ: Lawrence Erlbaum.
- Fang, Kai-Tai, Kotz, Samuel and Ng, Kai Wang (1990). *Symmetric Multivariate and Related Distributions*. London: Chapman and Hall.
- Fisher, Ronald A. (1922). The Goodness of Fit of Regression Formulae, and the Distribution of Regression Coefficients. *Journal of the Royal Statistical Society* 85(4), 597–612.
- Fisher, Ronald A. (1925). *Statistical Methods for Research Workers*. Edinburgh: Oliver and Boyd.
- Ghosh, B. K. (1973). Some Monotonicity Theorems for \( \chi^2 \), F and t Distributions with Applications. *Journal of the Royal Statistical Society. Series B* 35(3), 480–492.
- Hoenig, John M. and Heisey, Dennis M. (2001). The Abuse of Power: The Pervasive Fallacy of Power Calculations for Data Analysis. *The American Statistician* 55(1), 19–24.
- Lehmann, E. L. and Romano, Joseph P. (2005). *Testing Statistical Hypotheses*. 3rd edition. New York: Springer.
- Longley, James W. (1967). An Appraisal of Least Squares Programs for the Electronic Computer from the Point of View of the User. *Journal of the American Statistical Association* 62(319), 819–841.
- Neyman, Jerzy and Pearson, Egon S. (1933). On the Problem of the Most Efficient Tests of Statistical Hypotheses. *Philosophical Transactions of the Royal Society of London. Series A* 231, 289–337.
- Rao, C. Radhakrishna (1948). Large Sample Tests of Statistical Hypotheses Concerning Several Parameters with Applications to Problems of Estimation. *Mathematical Proceedings of the Cambridge Philosophical Society* 44(1), 50–57.
- Rencher, Alvin C. and Schaalje, G. Bruce (2008). *Linear Models in Statistics*. 2nd edition. Hoboken, NJ: Wiley.
- Scheffé, Henry (1959). *The Analysis of Variance*. New York: Wiley.
- Schuirmann, Donald J. (1987). A Comparison of the Two One-Sided Tests Procedure and the Power Approach for Assessing the Equivalence of Average Bioavailability. *Journal of Pharmacokinetics and Biopharmaceutics* 15(6), 657–680.
- Searle, Shayle R. (1971). *Linear Models*. New York: Wiley.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Sen, Ashish and Srivastava, Muni (1990). *Regression Analysis: Theory, Methods, and Applications*. New York: Springer.
- Snedecor, George W. (1934). *Calculation and Interpretation of Analysis of Variance and Covariance*. Ames, IA: Collegiate Press.
- Tang, P. C. (1938). The Power Function of the Analysis of Variance Tests with Tables and Illustrations of Their Use. *Statistical Research Memoirs* 2, 126–149.
- Wald, Abraham (1943). Tests of Statistical Hypotheses Concerning Several Parameters When the Number of Observations Is Large. *Transactions of the American Mathematical Society* 54(3), 426–482.
- Wasserstein, Ronald L. and Lazar, Nicole A. (2016). The ASA's Statement on p-Values: Context, Process, and Purpose. *The American Statistician* 70(2), 129–133.
- Wilks, Samuel S. (1938). The Large-Sample Distribution of the Likelihood Ratio for Testing Composite Hypotheses. *The Annals of Mathematical Statistics* 9(1), 60–62.
