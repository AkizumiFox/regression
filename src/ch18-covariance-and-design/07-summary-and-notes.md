# Summary and Notes

## Summary

::: {.idea}

1. Analysis of covariance is a projection problem. With design part \( \X \) (any rank) and covariates \( \Z \), the model space
   splits orthogonally as \( \C(\X)\dirsum\C((\I-\M)\Z) \). The covariate coefficients are estimated from within-group variation,
   \( \hat{\bgamma}=\mathbf{E}_{zz}^{-1}\mathbf{E}_{zy} \), and the design part is fitted to the adjusted response \( \y-\Z\hat{\bgamma} \). Tests compare
   residual sums of squares of nested models (@thm-dsn-ancova).

2. In a randomized experiment a covariate measured before randomization reduces the variance of treatment comparisons by the factor
   \( (1-\rho^2)\{1+1/(N-g-2)\} \) (@prp-dsn-precision), and it removes the conditional bias caused by chance imbalance
   ([Section 18.1](01-ancova.html)). A covariate affected by the treatment must not be adjusted for (same section).

3. Adjusted means are estimable functions \( \mu_k+\bz_0\T\bgamma \), estimated by \( \bar{y}_k-(\bar{\bz}_k-\bz_0)\T\hat{\bgamma} \). They are correlated,
   their variances grow with the distance of \( \bar{\bz}_k \) from \( \bz_0 \), and their differences do not depend on \( \bz_0 \) (@def-dsn-adjusted-means and @prp-dsn-adjusted-se). Scheffé's intervals are exact for all their contrasts; Tukey-type intervals are
   approximate.

4. The test of parallel slopes is a weighted test of homogeneity of the within-group slopes (@thm-dsn-slopes). When slopes differ, the
   treatment difference is a line in the covariate, to be reported with a band or a Johnson–Neyman region (@prp-dsn-johnson-neyman).

5. In a randomized complete block design, block effects cancel from every treatment contrast. The analysis is the additive two-way
   analysis (@thm-dsn-rcbd). Randomization alone makes the contrast estimates unbiased and the error mean square unbiased for their
   variance (@prp-dsn-randomization). The block experiment estimates the error variance a completely randomized design would have had (@prp-dsn-crd-efficiency).

6. A Latin square blocks in two crossed directions with mutually orthogonal row, column and treatment spaces. Its error consists of
   confounded interaction pieces (@thm-dsn-latin-square).

7. A balanced incomplete block design has information matrix \( (\lambda t/k)(\I-\mathbf{J}/t) \). Contrasts are estimated from the adjusted
   totals \( Q_i \) with variance \( (k/\lambda t)\sigma^2\sum c_i^2 \), and the efficiency factor relative to complete blocks is
   \( \lambda t/(rk) \) (@thm-dsn-bibd).

8. A lost observation is handled exactly by least squares on the remaining data. Filling it with its fitted value, found by Yates's formula,
   by iteration or by the covariate method, reproduces the exact estimates and residual sum of squares. The error degrees of freedom must be
   reduced, and tests must be computed by comparing models (@thm-dsn-missing and @cor-dsn-yates).

:::

## Notes and sources

**Analysis of covariance.**  Fisher introduced the analysis of covariance in the 1930s, in later editions of *Statistical Methods for
Research Workers* (Fisher 1925 and later editions), as a way of using a concomitant measurement to increase the precision of agricultural
experiments. Cochran's review (Cochran 1957) remains an excellent account of its uses, including the dangers of adjusting for a
covariate affected by treatment and of adjusting non-randomized groups. Rencher and Schaalje (2008, chapter 16) treat one-way, two-way and multi-covariate layouts and the test of homogeneous slopes. Christensen (2020, chapter 9) derives the analysis through the projection \( (\I-\M)\Z \), which is the route
taken in @thm-dsn-ancova, and uses the same device for missing data and for incomplete block designs. The case against
testing baseline balance in randomized trials was made forcefully by Senn (1994). The modern debate on regression adjustment in randomized
experiments, which asks whether adjustment can hurt when the model is wrong, was opened by Freedman (2008) and answered, for the
interacted estimator of @exr-dsn-interacted, by Lin (2013). [Chapter 25](../ch25-causal-interpretation/index.html) returns to it (@thm-cau-adjustment).

**Heterogeneous slopes.**  The region of significance for non-parallel lines is due to Johnson and Neyman (1936), who developed it for
comparisons of teaching methods with an ability covariate. Simultaneous versions use the Scheffé argument of
[Chapter 13](../ch13-multiplicity/index.html).

**Block designs and randomization.**  Randomization, blocking and the Latin square as a design are Fisher's, set out in *The Design of
Experiments* (Fisher 1935). The randomization moments of @prp-dsn-randomization and the view of the linear model as an approximation
justified by randomization are developed in Kempthorne (1952) and Hinkelmann and Kempthorne (2008). Cox (1958) remains the best short
account of the ideas behind designed experiments, and Bailey (2008) gives a modern treatment of block designs built on the same
orthogonal decompositions as this chapter. Relative efficiencies such as @eq-dsn-crd-variance are standard in the design literature, for instance Cochran and Cox
(1957). Christensen (2020, chapter 8) presents the randomized complete block and Latin square models as
additive multi-way layouts and discusses their assumptions.

**Incomplete blocks.**  Balanced incomplete block designs were introduced by Yates (1936), who also showed how to recover
inter-block information (Yates 1940). Fisher (1940) proved the inequality \( b\ge t \) of @lem-dsn-bibd. Christensen (2020, section 9.4)
derives the intra-block analysis as an analysis of covariance with the treatment indicators as covariates, as in the proof of @thm-dsn-bibd.

**Missing values.**  Yates (1933) derived the missing value formula @eq-dsn-yates and the rule of subtracting a degree of freedom for each
estimated value. Bartlett (1937) proposed the covariance method with indicator covariates, and Healy and Westmacott (1956) the iterative
method of @thm-dsn-missing(c), convenient for early computers. Dempster, Laird and Rubin (1977) placed the iteration within the EM
algorithm. Christensen (2020, section 9.3) treats missing observations through the covariance model.

## References

- Bailey, Rosemary A. (2008). *Design of Comparative Experiments*. Cambridge: Cambridge University Press.
- Bartlett, Maurice S. (1937). Some Examples of Statistical Methods of Research in Agriculture and Applied Biology. *Supplement to the Journal of the Royal Statistical Society* 4(2), 137–183.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Cochran, William G. (1957). Analysis of Covariance: Its Nature and Uses. *Biometrics* 13(3), 261–281.
- Cochran, William G. and Cox, Gertrude M. (1957). *Experimental Designs*. 2nd edition. New York: Wiley.
- Cox, David R. (1958). *Planning of Experiments*. New York: Wiley.
- Dempster, Arthur P., Laird, Nan M. and Rubin, Donald B. (1977). Maximum Likelihood from Incomplete Data via the EM Algorithm. *Journal of the Royal Statistical Society, Series B* 39(1), 1–38.
- Fisher, Ronald A. (1925). *Statistical Methods for Research Workers*. Edinburgh: Oliver and Boyd.
- Fisher, Ronald A. (1935). *The Design of Experiments*. Edinburgh: Oliver and Boyd.
- Fisher, Ronald A. (1940). An Examination of the Different Possible Solutions of a Problem in Incomplete Blocks. *Annals of Eugenics* 10(1), 52–75.
- Freedman, David A. (2008). On Regression Adjustments to Experimental Data. *Advances in Applied Mathematics* 40(2), 180–193.
- Healy, Michael and Westmacott, Michael (1956). Missing Values in Experiments Analysed on Automatic Computers. *Journal of the Royal Statistical Society, Series C (Applied Statistics)* 5(3), 203–206.
- Hinkelmann, Klaus and Kempthorne, Oscar (2008). *Design and Analysis of Experiments, Volume 1: Introduction to Experimental Design*. 2nd edition. Hoboken, NJ: Wiley.
- Johnson, Palmer O. and Neyman, Jerzy (1936). Tests of Certain Linear Hypotheses and Their Application to Some Educational Problems. *Statistical Research Memoirs* 1, 57–93.
- Kempthorne, Oscar (1952). *The Design and Analysis of Experiments*. New York: Wiley.
- Lin, Winston (2013). Agnostic Notes on Regression Adjustments to Experimental Data: Reexamining Freedman's Critique. *The Annals of Applied Statistics* 7(1), 295–318.
- Rencher, Alvin C. and Schaalje, G. Bruce (2008). *Linear Models in Statistics*. 2nd edition. Hoboken, NJ: Wiley.
- Senn, Stephen (1994). Testing for Baseline Balance in Clinical Trials. *Statistics in Medicine* 13(17), 1715–1726.
- Yates, Frank (1933). The Analysis of Replicated Experiments When the Field Results Are Incomplete. *Empire Journal of Experimental Agriculture* 1, 129–142.
- Yates, Frank (1936). Incomplete Randomized Blocks. *Annals of Eugenics* 7(2), 121–140.
- Yates, Frank (1940). The Recovery of Inter-Block Information in Balanced Incomplete Block Designs. *Annals of Eugenics* 10(4), 317–325.
