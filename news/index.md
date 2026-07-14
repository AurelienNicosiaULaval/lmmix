# Changelog

## lmmix 0.3.0

- Added fixed-band Toeplitz residual covariance through specifications
  such as `structure = "toep(4)"`, while retaining full `"toep"`
  behavior.
- Added an exact likelihood and covariance comparison with official PROC
  MIXED Example 79.6, including its line-source irrigation data and
  runnable SAS specification.
- Extended [`confint()`](https://rdrr.io/r/stats/confint.html) with
  `"beta_"` and `"theta_"` selectors. Covariance intervals use log and
  Fisher-z transformations to respect variance and correlation bounds.
- Added reproducible parametric-bootstrap likelihood-ratio tests for
  nested models through `anova(..., test = "parametric.bootstrap")`.
- Added a vignette on covariance intervals, boundary-aware comparisons,
  and inspection of the simulated null reference distribution.

## lmmix 0.2.1

- Added reproducible comparisons with official SAS PROC MIXED examples
  79.1, 79.2, and 79.5 from the SAS/STAT 14.3 documentation.
- Added the split-plot, repeated-measures growth, and
  random-coefficients data used by those examples.
- Added a runnable SAS validation program and machine-readable published
  targets for covariance parameters, likelihood criteria, fixed effects,
  BLUPs, and type III tests.
- Added package and repository citation metadata for academic reuse.
- Added CRAN submission preparation metadata and clarified the
  distinction between published SAS targets and a fresh SAS execution
  transcript.

## lmmix 0.2.0

- Added Kenward-Roger covariance adjustment, denominator degrees of
  freedom, and scaled multi-degree-of-freedom tests for REML fits.
- Added independent crossed or nested random-effect terms through a
  backward-compatible list-of-formulas interface.
- Added likelihood-ratio comparisons for nested models. Models that
  differ in fixed effects can be refitted automatically with ML, while
  covariance-model comparisons require ML and report the
  boundary-distribution caveat.
- Added deterministic optimizer fallback attempts and retained complete
  attempt diagnostics in fitted objects.
- Added `na.action` support for `na.omit`, `na.exclude`, and `na.fail`,
  including restoration of excluded rows in fitted values, residuals,
  and augmented data.
- New grouping levels now require `allow.new.levels = TRUE` for
  conditional prediction instead of silently receiving a zero random
  contribution.
- Added simultaneous Bonferroni confidence intervals for adjusted
  pairwise marginal-mean comparisons.
- Added independent validation against `lme4`, `lmerTest`, `mmrm`, and
  `emmeans` for the new 0.2.0 capabilities.

## lmmix 0.1.0

- Added explicit profiled ML and REML estimation for Gaussian linear
  mixed models.
- Added independent, compound-symmetric, AR(1), Toeplitz, and
  unstructured residual covariance matrices.
- Added random intercepts and slopes with an unstructured covariance
  matrix.
- Added Satterthwaite fixed-effect inference, type III tests, estimated
  marginal means, and pairwise contrasts.
- Added standard model, broom, and `emmeans` methods.
- Added
  [`plot.lmm()`](https://aureliennicosiaulaval.github.io/lmmix/reference/lmm-methods.md)
  diagnostics that return `ggplot2` objects for residual, normal Q-Q,
  and observed-versus-fitted displays.
- Added human-readable table printing without dots in displayed headings
  while preserving standard `broom` column names for programmatic use.
- Added the multilocation example from Milliken and Johnson (2009,
  Section 28.3) and validation tests against SAS PROC MIXED, `nlme`,
  `lmerTest`, and `mmrm`.
- Added primary theoretical references for REML, Satterthwaite
  inference, multi-degree-of-freedom tests, and estimated marginal
  means.
- Expanded the validation vignette with executable comparisons against
  `nlme`, `mmrm`, and `lmerTest`, plus every stored PROC MIXED target
  and its observed numerical difference.
- Reworked the README, vignettes, package overview, and method
  documentation to describe the implemented API, validation evidence,
  and current limits.
- Added combined-model convergence tests for every residual covariance
  structure.
