# Summary and Notes

## Summary

::: {.idea}

1. The error rates of a family, familywise, per-family and false discovery, are ordered (@def-mc-error-rates, @prp-mc-error-rate-order).
   Control should be strong; a gate such as the
   protected LSD gives only weak control (@def-mc-strong-control, @prp-mc-gatekeeper). Simultaneous
   intervals give strong control, including of sign errors (@prp-mc-intervals-to-tests).

2. Scheffé: a family of estimable functions closed under linear combination is a subspace of
   \( \C(\X) \), and its largest squared \( t \) statistic is \( q \) times an \( F(q,n-r) \) variable. The intervals
   with multiplier \( \sqrt{qF_\alpha(q,\nu)} \) are exact, even for functions chosen after seeing the data (@thm-mc-scheffe);
   the \( F \) test rejects iff some Scheffé interval excludes its hypothesized value (@cor-mc-scheffe-f);
   smaller subspaces give shorter intervals ([Section 13.2](02-scheffe.html)).

3. Tukey: in a balanced one-way layout the intervals \( \bar{y}_k-\bar{y}_l\pm q_\alpha(g,\nu)s/\sqrt m \) are exact
   for all pairs, extend to all contrasts, and are the shortest of their form (@def-mc-studentized-range, @thm-mc-tukey, @prp-mc-tukey-shortest).
   The Tukey–Kramer intervals @eq-mc-tukey-kramer are
   conservative.

4. Bonferroni needs no assumption on dependence and bounds the expected number of false rejections (@thm-mc-bonferroni);
   Šidák's level is valid for independent estimates, and a shared \( s \) only
   helps (@prp-mc-sidak). Closed testing gives strong control (@thm-mc-closed-testing), and Holm's
   step-down procedure, its Bonferroni case, rejects everything Bonferroni does and possibly more (@thm-mc-holm).

5. The Benjamini–Hochberg procedure controls the false discovery rate at \( \pi_0q \) for independent
   null \( p \)-values (@thm-mc-bh), and Benjamini and Yekutieli's results extend it to positively
   dependent (PRDS) and arbitrarily dependent \( p \)-values.

6. Fix the family first: planned lists call for Bonferroni, Holm or exact multivariate \( t \) values,
   all pairs for Tukey, data snooping for Scheffé, screening for the false discovery rate
   ([Section 13.6](06-choosing.html)).

:::

## Notes and sources

**Sources for this chapter.** Christensen (2020, chapter 5) treats multiple comparisons within the
projection framework used here, describing Scheffé's method through a test space and showing that
the \( F \) test rejects exactly when some single-degree-of-freedom hypothesis in that space does.
[Section 13.2](02-scheffe.html) follows that view. His chapter also compares the least significant difference,
Bonferroni, Tukey, Newman–Keuls and Duncan procedures and discusses the Fisherian and
Neyman–Pearson views of multiple testing. Agresti (2015, section 3.5) gives a concise account of
Bonferroni's and Tukey's methods, the Tukey–Kramer extension and the Benjamini–Hochberg procedure.
The examples, simulations and the treatment of error rates, Šidák's inequality and closed testing
are this book's own.

**Scheffé and Tukey.** Scheffé (1953) introduced the method for all contrasts, and his book
(Scheffé 1959) develops it for general linear hypotheses and relates it to the \( F \) test and the
confidence ellipsoid. Tukey's 1953 manuscript circulated for decades before it appeared in his
collected works (Tukey 1994). Kramer (1956) proposed the extension to unequal group sizes, and Hayter (1984)
proved Tukey's conjecture that these intervals are conservative. The case of estimates with general
correlation, where neither form is exact, is handled by computing multivariate \( t \) probabilities;
Hothorn, Bretz and Westfall (2008) describe a general framework and Bretz, Hothorn and Westfall
(2010) a practical guide. Dunnett (1955) gave the exact procedure for comparisons with a control.

**Bonferroni, Šidák and Holm.** The name commemorates Bonferroni's (1936) work on inequalities for
probabilities of unions. Dunn (1961) popularized its use
for simultaneous confidence intervals. Šidák (1967) proved the inequality for rectangular regions of
the multivariate normal distribution that underlies @prp-mc-sidak for correlated estimates. Holm
(1979) introduced the step-down procedure, and Marcus, Peritz and Gabriel (1976) the closure
principle, of which it is a case. Hochberg (1988) proposed the step-up version, whose validity rests
on the inequality of Simes (1986).

**False discovery rate.** Benjamini and Hochberg (1995) defined the false discovery rate and proved @thm-mc-bh.
The leave-one-out proof given here is a standard later simplification. Benjamini and Yekutieli (2001) proved control under positive regression dependence
and gave the version valid under arbitrary dependence. Storey (2002) introduced estimation of
\( \pi_0 \) and the \( q \)-value, and Benjamini, Krieger and Yekutieli (2006) gave adaptive procedures
with proven control.

**Books.** Miller (1981) is the classic account of simultaneous inference, and Hochberg and Tamhane
(1987) and Hsu (1996) are comprehensive monographs. All three treat the procedures of Duncan (1955)
and of Newman (1939) and Keuls (1952), and the closed step-down procedures that improve on Tukey's.
Gelman, Hill and Yajima (2012) give a Bayesian view of why hierarchical models reduce the need for
multiplicity adjustments.

## References

- Agresti, Alan (2015). *Foundations of Linear and Generalized Linear Models*. Hoboken, NJ: Wiley.
- Benjamini, Yoav and Hochberg, Yosef (1995). Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. *Journal of the Royal Statistical Society, Series B* 57(1), 289–300.
- Benjamini, Yoav, Krieger, Abba M. and Yekutieli, Daniel (2006). Adaptive Linear Step-up Procedures That Control the False Discovery Rate. *Biometrika* 93(3), 491–507.
- Benjamini, Yoav and Yekutieli, Daniel (2001). The Control of the False Discovery Rate in Multiple Testing under Dependency. *The Annals of Statistics* 29(4), 1165–1188.
- Bonferroni, Carlo E. (1936). Teoria statistica delle classi e calcolo delle probabilità. *Pubblicazioni del R. Istituto Superiore di Scienze Economiche e Commerciali di Firenze* 8, 3–62.
- Bretz, Frank, Hothorn, Torsten and Westfall, Peter (2010). *Multiple Comparisons Using R*. Boca Raton, FL: Chapman & Hall/CRC.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Duncan, David B. (1955). Multiple Range and Multiple F Tests. *Biometrics* 11(1), 1–42.
- Dunn, Olive Jean (1961). Multiple Comparisons among Means. *Journal of the American Statistical Association* 56(293), 52–64.
- Dunnett, Charles W. (1955). A Multiple Comparison Procedure for Comparing Several Treatments with a Control. *Journal of the American Statistical Association* 50(272), 1096–1121.
- Gelman, Andrew, Hill, Jennifer and Yajima, Masanao (2012). Why We (Usually) Don't Have to Worry about Multiple Comparisons. *Journal of Research on Educational Effectiveness* 5(2), 189–211.
- Hayter, Anthony J. (1984). A Proof of the Conjecture That the Tukey–Kramer Multiple Comparisons Procedure Is Conservative. *The Annals of Statistics* 12(1), 61–75.
- Hochberg, Yosef (1988). A Sharper Bonferroni Procedure for Multiple Tests of Significance. *Biometrika* 75(4), 800–802.
- Hochberg, Yosef and Tamhane, Ajit C. (1987). *Multiple Comparison Procedures*. New York: Wiley.
- Holm, Sture (1979). A Simple Sequentially Rejective Multiple Test Procedure. *Scandinavian Journal of Statistics* 6(2), 65–70.
- Hothorn, Torsten, Bretz, Frank and Westfall, Peter (2008). Simultaneous Inference in General Parametric Models. *Biometrical Journal* 50(3), 346–363.
- Hsu, Jason C. (1996). *Multiple Comparisons: Theory and Methods*. London: Chapman & Hall.
- Keuls, M. (1952). The Use of the “Studentized Range” in Connection with an Analysis of Variance. *Euphytica* 1, 112–122.
- Kramer, Clyde Young (1956). Extension of Multiple Range Tests to Group Means with Unequal Numbers of Replications. *Biometrics* 12(3), 307–310.
- Marcus, Ruth, Peritz, Eric and Gabriel, K. Ruben (1976). On Closed Testing Procedures with Special Reference to Ordered Analysis of Variance. *Biometrika* 63(3), 655–660.
- Miller, Rupert G. (1981). *Simultaneous Statistical Inference*. 2nd edition. New York: Springer.
- Newman, D. (1939). The Distribution of Range in Samples from a Normal Population, Expressed in Terms of an Independent Estimate of Standard Deviation. *Biometrika* 31(1/2), 20–30.
- Scheffé, Henry (1953). A Method for Judging All Contrasts in the Analysis of Variance. *Biometrika* 40(1/2), 87–104.
- Scheffé, Henry (1959). *The Analysis of Variance*. New York: Wiley.
- Šidák, Zbyněk (1967). Rectangular Confidence Regions for the Means of Multivariate Normal Distributions. *Journal of the American Statistical Association* 62(318), 626–633.
- Simes, R. John (1986). An Improved Bonferroni Procedure for Multiple Tests of Significance. *Biometrika* 73(3), 751–754.
- Storey, John D. (2002). A Direct Approach to False Discovery Rates. *Journal of the Royal Statistical Society, Series B* 64(3), 479–498.
- Tukey, John W. (1994). The Problem of Multiple Comparisons (1953 manuscript). In H. I. Braun (ed.), *The Collected Works of John W. Tukey, Volume VIII: Multiple Comparisons, 1948–1983*. New York: Chapman & Hall.
