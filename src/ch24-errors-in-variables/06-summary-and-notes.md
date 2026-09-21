# Summary and Notes

## Summary

::: {.idea}

1. Under classical, nondifferential error, least squares estimates the best linear predictor from the *reading*.
   A simple-regression slope is attenuated by the reliability \( \lambda \), whatever \( n \) (@thm-eiv-attenuation).

2. Exactly measured covariates absorb part of the effect and can look important when they are not (@exm-eiv-leak).
   Error in the response only adds variance (@prp-eiv-response). The test of no effect of the mismeasured
   regressor keeps its level; tests on correlated, correctly measured covariates do not (@prp-eiv-tests).

3. Attenuation is the same with fixed true values (@prp-eiv-functional). The normal structural model identifies the
   slope only up to the interval between the direct and reverse regressions (@prp-eiv-nonidentified). Berkson error
   causes no bias in a linear model (@prp-eiv-berkson).

4. A known or replicate-estimated error covariance gives the moment correction, with a sandwich covariance (@thm-eiv-correction);
   a known variance ratio gives Deming regression (@prp-eiv-deming).

5. IV and two-stage least squares are consistent under exogeneity and relevance (@lem-eiv-2sls, @thm-eiv-iv); with weak
   instruments they have no mean, drift towards least squares, and Wald intervals undercover (@prp-eiv-weak).

6. Linear regression calibration is exact and equals the moment correction (@prp-eiv-calibration). SIMEX extrapolates
   the naive estimate to zero error; in the linear model the rational form matches the limit exactly, so the rational extrapolant is consistent under
   regularity conditions for its fit (@prp-eiv-simex).

:::

## Notes and sources

**Sources.** None of the six source books treats errors in variables at length. Seber and Lee (2003, section 9.6)
surveys random regressors, fixed regressors measured with error, the structural identification problem, rounding
errors and Berkson's model, covered here in @thm-eiv-attenuation, @prp-eiv-functional, @prp-eiv-nonidentified,
@exr-eiv-rounding and @prp-eiv-berkson. Fuller (1987) is the classical treatment of linear measurement error models,
including moment estimators, their asymptotics and small-sample modifications. Carroll, Ruppert, Stefanski and
Crainiceanu (2006) is the standard reference for nonlinear models, regression calibration and SIMEX, and discusses
the choice of SIMEX extrapolant. All data in the chapter are simulated.

**Attenuation and identification.** Spearman (1904) introduced the correction for attenuation in mental testing.
Reiersøl (1950) proved that nonnormality of the true regressor restores identification. The functional model's
incidental parameters are an instance of Neyman and Scott (1948). Gleser and Hwang (1987) proved that
errors-in-variables confidence sets with guaranteed coverage have infinite expected length. Berkson (1950)
distinguished controlled regressors from classical error.

**Orthogonal and Deming regression.** Orthogonal regression goes back to Adcock in the 1870s; the version with a
known variance ratio carries Deming's name from Deming (1943). In numerical analysis it is total least squares
(Golub and Van Loan, 1980). Wald (1940) proposed the grouping estimator of @exr-eiv-wald-grouping.

**Instrumental variables.** Instrumental variables go back to Wright (1928). Two-stage least squares was developed in
the 1950s by Theil and by Basmann (1957), among others, and Sargan (1958) gave the general theory, including the test
of overidentifying restrictions. Anderson and Rubin (1949) introduced the test of @exr-eiv-anderson-rubin. Bound,
Jaeger and Baker (1995) showed the practical danger of weak instruments, and Staiger and Stock (1997) gave the
asymptotics behind the first-stage \( F \) rule of thumb. Angrist and Pischke (2009) give a modern treatment, including the causal use of IV taken up in [Chapter 25](../ch25-causal-interpretation/index.html).

**Regression calibration and SIMEX.** Regression calibration developed in epidemiology in the 1980s to correct
relative risks for exposure measurement error. SIMEX is due to Cook and Stefanski (1994).

## References

- Anderson, T. W. and Rubin, Herman (1949). Estimation of the Parameters of a Single Equation in a Complete System of Stochastic Equations. *The Annals of Mathematical Statistics* 20(1), 46–63.
- Angrist, Joshua D. and Pischke, Jörn-Steffen (2009). *Mostly Harmless Econometrics: An Empiricist's Companion*. Princeton, NJ: Princeton University Press.
- Basmann, R. L. (1957). A Generalized Classical Method of Linear Estimation of Coefficients in a Structural Equation. *Econometrica* 25(1), 77–83.
- Berkson, Joseph (1950). Are There Two Regressions? *Journal of the American Statistical Association* 45(250), 164–180.
- Bound, John, Jaeger, David A. and Baker, Regina M. (1995). Problems with Instrumental Variables Estimation When the Correlation Between the Instruments and the Endogenous Explanatory Variable Is Weak. *Journal of the American Statistical Association* 90(430), 443–450.
- Carroll, Raymond J., Ruppert, David, Stefanski, Leonard A. and Crainiceanu, Ciprian M. (2006). *Measurement Error in Nonlinear Models: A Modern Perspective*. 2nd edition. Boca Raton, FL: Chapman & Hall/CRC.
- Cook, J. R. and Stefanski, L. A. (1994). Simulation-Extrapolation Estimation in Parametric Measurement Error Models. *Journal of the American Statistical Association* 89(428), 1314–1328.
- Deming, W. Edwards (1943). *Statistical Adjustment of Data*. New York: Wiley.
- Fuller, Wayne A. (1987). *Measurement Error Models*. New York: Wiley.
- Gleser, Leon Jay and Hwang, Jiunn Tzon (1987). The Nonexistence of 100(1−α)% Confidence Sets of Finite Expected Diameter in Errors-in-Variables and Related Models. *The Annals of Statistics* 15(4), 1351–1362.
- Golub, Gene H. and Van Loan, Charles F. (1980). An Analysis of the Total Least Squares Problem. *SIAM Journal on Numerical Analysis* 17(6), 883–893.
- Neyman, J. and Scott, Elizabeth L. (1948). Consistent Estimates Based on Partially Consistent Observations. *Econometrica* 16(1), 1–32.
- Reiersøl, Olav (1950). Identifiability of a Linear Relation between Variables Which Are Subject to Error. *Econometrica* 18(4), 375–389.
- Sargan, J. D. (1958). The Estimation of Economic Relationships Using Instrumental Variables. *Econometrica* 26(3), 393–415.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Spearman, C. (1904). The Proof and Measurement of Association between Two Things. *The American Journal of Psychology* 15(1), 72–101.
- Staiger, Douglas and Stock, James H. (1997). Instrumental Variables Regression with Weak Instruments. *Econometrica* 65(3), 557–586.
- Wald, Abraham (1940). The Fitting of Straight Lines if Both Variables Are Subject to Error. *The Annals of Mathematical Statistics* 11(3), 284–300.
- Wright, Philip G. (1928). *The Tariff on Animal and Vegetable Oils*. New York: Macmillan.
