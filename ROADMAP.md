# lmmix roadmap

This document centralizes the current scope, known limitations, validation
boundaries, and development priorities for `lmmix`. It describes version
0.3.0 as released on 2026-07-14.

The roadmap is informative rather than a promise of release dates. A feature
is considered implemented only when it appears in `NEWS.md`, is documented,
and is covered by automated tests.

## Release history

The entries below summarize the role of each released or documented version.
The complete itemized changelog is available in
[NEWS](https://aureliennicosiaulaval.github.io/lmmix/news/index.html).

### 0.1.0: core modeling and inference

Version 0.1.0 established the initial package implementation on 2026-07-13.
It introduced explicit profiled ML and REML estimation, the five original
residual covariance structures, random intercepts and slopes, Satterthwaite
inference, type III tests, estimated marginal means, standard model and
`broom` methods, `emmeans` interoperability, and `ggplot2` diagnostics. It also
established the theoretical and numerical validation documentation.

This initial development version was not published as a GitHub release and
does not have a corresponding Git tag. Its history begins with the
[initial implementation](https://github.com/AurelienNicosiaULaval/lmmix/commit/73e8b163241e42e713bf1e519b824fa9ab16fdeb).

### 0.2.0: broader inference and model handling

[Version 0.2.0](https://github.com/AurelienNicosiaULaval/lmmix/releases/tag/v0.2.0),
published on 2026-07-14, added Kenward-Roger inference, multiple independent
random-effect terms, nested-model likelihood-ratio comparisons, deterministic
optimizer fallbacks, explicit missing-data policies, guarded prediction for
new grouping levels, simultaneous Bonferroni intervals, and independent
validation of these capabilities.

### 0.2.1: external validation and release metadata

[Version 0.2.1](https://github.com/AurelienNicosiaULaval/lmmix/releases/tag/v0.2.1),
published on 2026-07-14, consolidated validation against official PROC MIXED
examples 79.1, 79.2, and 79.5. It added the corresponding data, runnable SAS
program, machine-readable published targets, academic citation metadata, and
release-preparation documentation.

### 0.3.0: boundary-aware comparisons

[Version 0.3.0](https://github.com/AurelienNicosiaULaval/lmmix/releases/tag/v0.3.0),
published on 2026-07-14, added fixed-band Toeplitz covariance, validation
against PROC MIXED example 79.6, transformed covariance intervals, and
parametric-bootstrap likelihood-ratio tests for boundary-sensitive model
comparisons.

## Current production scope

Version 0.3.0 fits univariate Gaussian linear mixed models by explicitly
optimizing the ML or REML criterion. It supports:

- one or more independent random-effect terms, including correlated random
  slopes within a term;
- independent, compound-symmetric, AR(1), full Toeplitz, fixed-band
  `toep(k)`, and unstructured residual covariance;
- combined random effects and correlated residuals;
- Satterthwaite, Kenward-Roger, and residual denominator degrees of freedom;
- type III tests, estimated marginal means, contrasts, BLUPs, prediction,
  covariance intervals, and nested-model comparisons;
- parametric-bootstrap likelihood-ratio tests for boundary-sensitive model
  comparisons.

## Current limitations

### Model scope

- Responses must be univariate, continuous, and Gaussian.
- Generalized responses, multivariate responses, and nonlinear mixed models
  are not implemented.
- One repeated-measures structure can be specified per fitted model.
- Externally fixed multivariate G and R matrices are not accepted through the
  public API.
- Missing responses and covariates can be omitted or rejected, but missing
  covariates are not imputed.

### Computation and scale

- The marginal covariance matrix is assembled and factorized as a dense
  matrix at every objective evaluation.
- The package is intended for small and moderate data sets.
- Large-scale sparse algorithms, distributed fitting, and automatic
  parallelization are not implemented.
- Parametric-bootstrap comparisons refit both models for every simulation and
  can therefore be computationally expensive.

### Inference

- Containment denominator degrees of freedom are not implemented.
- Covariance confidence intervals are local Wald intervals. Variances use a
  log transformation and correlations use a Fisher-z transformation. They are
  not profile-likelihood intervals.
- Parametric-bootstrap p-values have Monte Carlo error. The installed vignette
  uses 49 simulations only as a fast demonstration; substantive analyses
  should use a larger `nsim`, commonly at least 999.
- Deletion and influence diagnostics comparable to PROC MIXED examples 79.7
  and 79.8 are not implemented.
- Conditional predictions for unseen grouping levels require
  `allow.new.levels = TRUE`; their random-effect contribution is then zero.

## Validation boundaries

- Automated comparisons cover overlapping models from `nlme`, `lme4`,
  `lmerTest`, `mmrm`, and `emmeans`.
- Published PROC MIXED targets cover official examples 79.1, 79.2, 79.5, and
  79.6, plus the stored multilocation comparison.
- The official SAS values are transcribed from published output tables. The
  repository contains a runnable SAS program but no archived log from a fresh
  SAS execution.
- The combined random-effect and correlated-residual Kenward-Roger case has
  structural tests but no independent SAS execution target.
- Generalized responses and large-scale performance are not validated because
  they are outside the current model scope.

Detailed numerical evidence is available in `VALIDATION.md` and in the
[validation vignette](https://aureliennicosiaulaval.github.io/lmmix/articles/validation.html).

## Development priorities

### Validation and reproducibility

- Archive a reproducible SAS log and output for the supported PROC MIXED
  examples.
- Add an independent target for the combined-model Kenward-Roger calculation.
- Add performance benchmarks with clearly reported sample sizes and covariance
  structures.

### Inference

- Evaluate profile-likelihood or bootstrap confidence intervals for covariance
  parameters.
- Evaluate containment degrees of freedom for designs where that convention is
  required for direct PROC MIXED reproduction.
- Add influence and deletion diagnostics with documented statistical
  definitions.
- Evaluate parallel execution for parametric-bootstrap comparisons while
  preserving reproducible random-number streams.

### Computation

- Exploit repeated-measures blocks more directly during covariance assembly and
  factorization.
- Evaluate sparse or blockwise algorithms before extending the package to
  larger data sets.
- Add explicit performance regression tests before making scalability claims.

## Items not currently scheduled

Generalized linear mixed models, nonlinear mixed models, and multivariate
responses would require substantial changes to the likelihood and inference
layers. They are not assigned to a planned release at this time.

## Tracking changes

- Released changes are recorded in the
  [NEWS changelog](https://aureliennicosiaulaval.github.io/lmmix/news/index.html).
- Numerical evidence and tolerances are recorded in `VALIDATION.md`.
- Proposed work can be discussed through
  [GitHub issues](https://github.com/AurelienNicosiaULaval/lmmix/issues).
