# lmmix roadmap

This document centralizes the current scope, known limitations,
validation boundaries, and development priorities for `lmmix`. It
describes version 0.4.0 as released on 2026-07-15.

The roadmap is informative rather than a promise of release dates. A
feature is considered implemented only when it appears in `NEWS.md`, is
documented, and is covered by automated tests.

## Release history

The entries below summarize the role of each released or documented
version. The complete itemized changelog is available in
[NEWS](https://aureliennicosiaulaval.github.io/lmmix/news/index.html).

### 0.1.0: core modeling and inference

Version 0.1.0 established the initial package implementation on
2026-07-13. It introduced explicit profiled ML and REML estimation, the
five original residual covariance structures, random intercepts and
slopes, Satterthwaite inference, type III tests, estimated marginal
means, standard model and `broom` methods, `emmeans` interoperability,
and `ggplot2` diagnostics. It also established the theoretical and
numerical validation documentation.

This initial development version was not published as a GitHub release
and does not have a corresponding Git tag. Its history begins with the
[initial
implementation](https://github.com/AurelienNicosiaULaval/lmmix/commit/73e8b163241e42e713bf1e519b824fa9ab16fdeb).

### 0.2.0: broader inference and model handling

[Version
0.2.0](https://github.com/AurelienNicosiaULaval/lmmix/releases/tag/v0.2.0),
published on 2026-07-14, added Kenward-Roger inference, multiple
independent random-effect terms, nested-model likelihood-ratio
comparisons, deterministic optimizer fallbacks, explicit missing-data
policies, guarded prediction for new grouping levels, simultaneous
Bonferroni intervals, and independent validation of these capabilities.

### 0.2.1: external validation and release metadata

[Version
0.2.1](https://github.com/AurelienNicosiaULaval/lmmix/releases/tag/v0.2.1),
published on 2026-07-14, consolidated validation against official PROC
MIXED examples 79.1, 79.2, and 79.5. It added the corresponding data,
runnable SAS program, machine-readable published targets, academic
citation metadata, and release-preparation documentation.

### 0.3.0: boundary-aware comparisons

[Version
0.3.0](https://github.com/AurelienNicosiaULaval/lmmix/releases/tag/v0.3.0),
published on 2026-07-14, added fixed-band Toeplitz covariance,
validation against PROC MIXED example 79.6, transformed covariance
intervals, and parametric-bootstrap likelihood-ratio tests for
boundary-sensitive model comparisons.

### 0.4.0: blockwise likelihood and standard model extensions

[Version
0.4.0](https://github.com/AurelienNicosiaULaval/lmmix/releases/tag/v0.4.0),
published on 2026-07-15, introduced connected-component block
factorization for the likelihood, cached Satterthwaite derivatives,
simulation, model updating, prediction intervals, partial random-effect
prediction, known precision weights, offsets, and custom contrasts. It
also standardized the
[`ranef()`](https://aureliennicosiaulaval.github.io/lmmix/reference/ranef.md)
return type and normal [`print()`](https://rdrr.io/r/base/print.html)
output.

## Current production scope

Version 0.4.0 fits univariate Gaussian linear mixed models by explicitly
optimizing the ML or REML criterion. It supports:

- one or more independent random-effect terms, including correlated
  random slopes within a term;
- independent, compound-symmetric, AR(1), full Toeplitz, fixed-band
  `toep(k)`, and unstructured residual covariance;
- combined random effects and correlated residuals;
- Satterthwaite, Kenward-Roger, and residual denominator degrees of
  freedom;
- type III tests, estimated marginal means, contrasts, BLUPs,
  prediction, covariance intervals, and nested-model comparisons;
- parametric-bootstrap likelihood-ratio tests for boundary-sensitive
  model comparisons;
- marginal simulation, model updating, prediction intervals, known
  relative precision weights, offsets, and custom fixed-effect
  contrasts.

## Current limitations

### Model scope

- Responses must be univariate, continuous, and Gaussian.
- Generalized responses, multivariate responses, and nonlinear mixed
  models are not implemented.
- One repeated-measures structure can be specified per fitted model.
- Externally fixed multivariate G and R matrices are not accepted
  through the public API.
- Missing responses and covariates can be omitted or rejected, but
  missing covariates are not imputed.
- Rank-deficient fixed-effect model matrices are rejected. Automatically
  dropping aliased columns is deferred until estimability and type III
  contrast behavior can be defined and tested together.
- Residual variance strata estimated from the data, analogous to
  [`nlme::varIdent()`](https://rdrr.io/pkg/nlme/man/varIdent.html), are
  not implemented. `weights` represent known relative precisions and are
  not estimated.

### Computation and scale

- Independent connected covariance components are assembled and
  factorized separately during likelihood evaluation.
- The fitted object still stores a complete dense marginal covariance
  matrix, and some inference calculations remain dense.
- Fully crossed designs can form one large connected component.
  Practical cost depends primarily on the largest component, the
  covariance parameter count, and the requested inference.
- Sparse Woodbury algorithms, distributed fitting, and automatic
  parallelization are not implemented. No fixed maximum sample size is
  claimed.
- Parametric-bootstrap comparisons refit both models for every
  simulation and can therefore be computationally expensive.

### Inference

- Containment denominator degrees of freedom are not implemented.
- Covariance confidence intervals are local Wald intervals. Variances
  use a log transformation and correlations use a Fisher-z
  transformation. They are not profile-likelihood intervals.
- Parametric-bootstrap p-values have Monte Carlo error. The installed
  vignette uses 49 simulations only as a fast demonstration; substantive
  analyses should use a larger `nsim`, commonly at least 999.
- Deletion and influence diagnostics comparable to PROC MIXED examples
  79.7 and 79.8 are not implemented.
- Conditional predictions for unseen grouping levels require
  `allow.new.levels = TRUE`; their random-effect contribution is then
  zero.
- Prediction confidence intervals condition on included empirical random
  effects. Prediction intervals add residual variance but not BLUP
  estimation uncertainty.
- The Kenward-Roger implementation uses finite-difference covariance
  derivatives. Analytic derivatives are not yet implemented.

## Validation boundaries

- Automated comparisons cover overlapping models from `nlme`, `lme4`,
  `lmerTest`, `mmrm`, and `emmeans`.
- Published PROC MIXED targets cover official examples 79.1, 79.2, 79.5,
  and 79.6, plus the stored multilocation comparison.
- The official SAS values are transcribed from published output tables.
  The repository contains a runnable SAS program but no archived log
  from a fresh SAS execution.
- The combined random-effect and correlated-residual Kenward-Roger case
  has structural tests but no independent SAS execution target.
- Generalized responses and large-scale performance are not validated
  because they are outside the current model scope.

Detailed numerical evidence is available in `VALIDATION.md` and in the
[validation
vignette](https://aureliennicosiaulaval.github.io/lmmix/articles/validation.html).

## Development priorities

### Validation and reproducibility

- Archive a reproducible SAS log and output for the supported PROC MIXED
  examples.
- Add an independent target for the combined-model Kenward-Roger
  calculation.
- Add performance benchmarks with clearly reported sample sizes and
  covariance structures.

### Inference

- Evaluate profile-likelihood or bootstrap confidence intervals for
  covariance parameters.
- Evaluate containment degrees of freedom for designs where that
  convention is required for direct PROC MIXED reproduction.
- Add influence and deletion diagnostics with documented statistical
  definitions.
- Evaluate parallel execution for parametric-bootstrap comparisons while
  preserving reproducible random-number streams.
- Implement estimability-aware handling of aliased fixed-effect columns
  and define its consequences for type III tests and marginal means.
- Evaluate estimated residual variance strata with an explicit
  covariance parameterization and validation against `nlme` or PROC
  MIXED.

### Computation

- Evaluate a sparse Woodbury formulation for large crossed random-effect
  designs.
- Derive and test analytic covariance derivatives for supported residual
  and random-effect structures, including their use in Kenward-Roger
  inference.
- Add explicit performance regression tests before making scalability
  claims.

## Items not currently scheduled

Generalized linear mixed models, nonlinear mixed models, and
multivariate responses would require substantial changes to the
likelihood and inference layers. They are not assigned to a planned
release at this time.

## Tracking changes

- Released changes are recorded in the [NEWS
  changelog](https://aureliennicosiaulaval.github.io/lmmix/news/index.html).
- Numerical evidence and tolerances are recorded in `VALIDATION.md`.
- Proposed work can be discussed through [GitHub
  issues](https://github.com/AurelienNicosiaULaval/lmmix/issues).
