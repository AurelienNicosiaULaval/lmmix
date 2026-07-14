# lmmix roadmap

This document centralizes the current scope, known limitations, validation
boundaries, and development priorities for `lmmix`. It describes version
0.3.0 as released on 2026-07-14.

The roadmap is informative rather than a promise of release dates. A feature
is considered implemented only when it appears in `NEWS.md`, is documented,
and is covered by automated tests.

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

- Released changes are recorded in `NEWS.md`.
- Numerical evidence and tolerances are recorded in `VALIDATION.md`.
- Proposed work can be discussed through
  [GitHub issues](https://github.com/AurelienNicosiaULaval/lmmix/issues).
