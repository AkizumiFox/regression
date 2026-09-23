# Summary and Notes

## Summary

::: {.idea}

1. Confidence sets come from pivots, and are the sets of values that a family of level-\( \alpha \) tests accepts
   (@def-ci-confidence-set, @prp-ci-duality).

2. For estimable \( \blambda\T\bbeta \), \( \blambda\T\hbeta\pm t_{n-r,\alpha/2}\,s\sqrt{\blambda\T\G\blambda} \) is exact and free of the
   generalized inverse (@thm-ci-estimable-interval). A non-estimable function has no confidence set of finite expected length
   (@prp-ci-nonestimable). The interval for \( \sigma^2 \) is exact under normality and fragile without it (@prp-ci-sigma).

3. The confidence ellipsoid for \( \bLambda\T\bbeta \) is the acceptance region of the \( F \) test, shaped by the covariance of the estimates
   (@thm-ci-ellipsoid); all such ellipsoids are images of one ball around \( \hY \) (@prp-ci-mean-ball). Its shadows are Scheffé-type intervals
   (@lem-ci-cauchy-schwarz, @cor-ci-shadows).

4. The prediction interval \( \hat{Y}_0\pm t_{n-r,\alpha/2}\,s\sqrt{1+h_0} \) never shrinks below the noise (@thm-ci-prediction-interval). The
   leverage \( h_0 \) is a Mahalanobis distance that reveals hidden extrapolation (@prp-ci-new-leverage). Without normality prediction
   intervals stay wrong as \( n\to\infty \) (@prp-ci-pi-limit).

5. The Working–Hotelling band, with multiplier \( \sqrt{kF_\alpha(k,n-r)} \), covers the regression surface over a \( k \)-dimensional family
   exactly; \( k=2 \) for a straight line (@thm-ci-working-hotelling).

6. Fieller sets for ratios, and calibration sets, are exact, and are intervals, two rays or the whole line (@prp-ci-fieller,
   @thm-ci-calibration). The classical calibration estimator has no mean (@prp-ci-ratio-no-mean), and unbounded sets cannot be avoided
   (@prp-ci-calibration-unbounded).

7. Conjugate credible sets and predictive intervals have the frequentist shapes (@lem-ci-mvt, @thm-ci-bayes-credible) and equal them
   under the flat prior (@cor-ci-flat-agreement). Their coverage is nominal on average over the prior, not pointwise
   (@prp-ci-average-coverage).

:::

## Notes and sources

**Confidence intervals and duality.** Confidence intervals, and their construction by inverting tests, are due to Neyman (1937). Lehmann and
Romano (2005) develop the duality and its optimality consequences. The two one-sided tests
procedure of @exr-ci-equivalence is due to Schuirmann (1987). The shortest interval for a normal variance (@exr-ci-shortest-sigma) was
worked out by Tate and Klett (1959).

**The source books.** The chapter covers Seber and Lee (2003, chapter 5 and sections 6.1–6.5), Rencher and Schaalje
(2008, sections 8.6–8.7) and Sen and Srivastava (1990, section 3.8). These books state the intervals for full-rank models; here they
are stated for estimable functions in models of any rank, as in Christensen (2020). Rencher and Schaalje (2008, section 8.7) derive the \( F \) test as a likelihood ratio test, which
[Chapter 11](../ch11-general-linear-hypothesis/index.html) covers (@thm-glh-lrt). Hahn and Meeker (1991) is a practical guide to
confidence, prediction and tolerance intervals.

**Bands.** The band for a straight line is due to Working and Hotelling (1929), and the general projection argument to
Scheffé (1953). Exact bands over finite intervals are due to Wynn and Bloomfield (1971) and constant-width bands to Gafarian (1964).
Miller (1981) surveys the early literature and Liu (2010) gives a modern book-length treatment.

**Normality and prediction.** The warning about
conditional coverage after a preliminary test is from Olshen (1973).

**Ratios and calibration.** Fieller's method came from bioassay and is set out in general form in Fieller (1954). Krutchkoff (1967)
proposed the inverse estimator and Hoadley (1970) gave its Bayesian interpretation. Osborne (1991) reviews calibration, and Brown (1993)
treats the multivariate case. Gleser and Hwang (1987) proved the general form of @prp-ci-calibration-unbounded.

**Bayesian intervals.** Box and Tiao (1973) and Gelman et al. (2013) treat credible and predictive intervals for the normal
linear model under conjugate and flat priors.

**Data.** The murder-rate data are the 2009 US state data of statsmodels; the stack loss data are from Brownlee (1965); the
assay is simulated.

## References

- Box, George E. P. and Tiao, George C. (1973). *Bayesian Inference in Statistical Analysis*. Reading, MA: Addison-Wesley.
- Brown, Philip J. (1993). *Measurement, Regression, and Calibration*. Oxford: Clarendon Press.
- Brownlee, Kenneth A. (1965). *Statistical Theory and Methodology in Science and Engineering*. 2nd edition. New York: Wiley.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Fieller, Edgar C. (1954). Some Problems in Interval Estimation. *Journal of the Royal Statistical Society, Series B* 16(2), 175–185.
- Gafarian, A. V. (1964). Confidence Bands in Straight Line Regression. *Journal of the American Statistical Association* 59(305), 182–213.
- Gelman, Andrew, Carlin, John B., Stern, Hal S., Dunson, David B., Vehtari, Aki and Rubin, Donald B. (2013). *Bayesian Data Analysis*. 3rd edition. Boca Raton, FL: CRC Press.
- Gleser, Leon Jay and Hwang, Jiunn Tzon (1987). The Nonexistence of \( 100(1-\alpha) \)% Confidence Sets of Finite Expected Diameter in Errors-in-Variables and Related Models. *The Annals of Statistics* 15(4), 1351–1362.
- Hahn, Gerald J. and Meeker, William Q. (1991). *Statistical Intervals: A Guide for Practitioners*. New York: Wiley.
- Hoadley, Bruce (1970). A Bayesian Look at Inverse Linear Regression. *Journal of the American Statistical Association* 65(329), 356–369.
- Krutchkoff, Richard G. (1967). Classical and Inverse Regression Methods of Calibration. *Technometrics* 9(3), 425–439.
- Lehmann, Erich L. and Romano, Joseph P. (2005). *Testing Statistical Hypotheses*. 3rd edition. New York: Springer.
- Liu, Wei (2010). *Simultaneous Inference in Regression*. Boca Raton, FL: Chapman & Hall/CRC.
- Miller, Rupert G. (1981). *Simultaneous Statistical Inference*. 2nd edition. New York: Springer.
- Neyman, Jerzy (1937). Outline of a Theory of Statistical Estimation Based on the Classical Theory of Probability. *Philosophical Transactions of the Royal Society of London. Series A* 236(767), 333–380.
- Olshen, Richard A. (1973). The Conditional Level of the F-Test. *Journal of the American Statistical Association* 68(343), 692–698.
- Osborne, Christine (1991). Statistical Calibration: A Review. *International Statistical Review* 59(3), 309–336.
- Rencher, Alvin C. and Schaalje, G. Bruce (2008). *Linear Models in Statistics*. 2nd edition. Hoboken, NJ: Wiley.
- Scheffé, Henry (1953). A Method for Judging All Contrasts in the Analysis of Variance. *Biometrika* 40(1/2), 87–104.
- Schuirmann, Donald J. (1987). A Comparison of the Two One-Sided Tests Procedure and the Power Approach for Assessing the Equivalence of Average Bioavailability. *Journal of Pharmacokinetics and Biopharmaceutics* 15(6), 657–680.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Sen, Ashish and Srivastava, Muni (1990). *Regression Analysis: Theory, Methods, and Applications*. New York: Springer.
- Tate, Robert F. and Klett, G. W. (1959). Optimal Confidence Intervals for the Variance of a Normal Distribution. *Journal of the American Statistical Association* 54(287), 674–682.
- Working, Holbrook and Hotelling, Harold (1929). Applications of the Theory of Error to the Interpretation of Trends. *Journal of the American Statistical Association* 24(165A), 73–85.
- Wynn, Henry P. and Bloomfield, Peter (1971). Simultaneous Confidence Bands in Regression Analysis. *Journal of the Royal Statistical Society, Series B* 33(2), 202–217.
