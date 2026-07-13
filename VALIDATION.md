# Validation scope

## Evidence layers

The package uses three complementary validation layers:

1.  structural tests for covariance parameterization, positive
    definiteness, dimensions, missing-data handling, input validation,
    methods, and output classes;
2.  numerical comparisons with independent R implementations in
    overlapping model classes;
3.  stored numerical regression targets for a combined random-effect and
    AR(1) residual model specified in PROC MIXED.

No single external R package covers every covariance and inference
component implemented by `lmmix`. Claims are therefore restricted to the
quantities tested in each overlapping case.

## Independent R comparisons

| Model | Reference | Quantities checked | Tolerance |
|:---|:---|:---|:---|
| Random-intercept REML | [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) | Fixed effects and log-likelihood | absolute `1e-6` |
| Random-intercept REML | [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) | Covariance parameters | absolute `1e-4` |
| Random-slope REML | [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) | Fixed effects and log-likelihood | absolute `1e-5` |
| Random-slope REML | [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) | Random covariance and residual variance | absolute `1e-4` |
| Random-intercept ML | [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) | Fixed effects and log-likelihood | `1e-6` |
| Marginal CS and AR(1) | [`nlme::gls()`](https://rdrr.io/pkg/nlme/man/gls.html) | Fixed effects, residual parameters, log-likelihood | `1e-5` |
| Marginal AR(1) | [`mmrm::mmrm()`](https://openpharma.github.io/mmrm/latest-tag/reference/mmrm.html) | Fixed effects and log-likelihood | absolute `1e-5` |
| Marginal AR(1) | [`mmrm::mmrm()`](https://openpharma.github.io/mmrm/latest-tag/reference/mmrm.html) | Residual covariance matrix | absolute `1e-4` |
| Satterthwaite coefficients | [`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html) | Estimates and standard errors | absolute `1e-5` |
| Satterthwaite coefficients | [`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html) | Denominator df | absolute `5e-3` |
| Type III Satterthwaite tests | [`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html) | Numerator df, denominator df, F statistics | absolute `1e-6`, `5e-3`, and `5e-3` |

The combined-model tests fit `id`, `cs`, `ar1`, `toep`, and `un`
residual structures together with a random center intercept. Every fit
must converge with a positive-definite likelihood Hessian and marginal
covariance matrix.

## Multilocation data

The included `multicentre` data reproduce the multilocation
repeated-measures example in Section 28.3 of Milliken and Johnson
(2009). The data contain 153 rows, including 28 missing responses, so
the model below uses 125 complete observations.

Milliken, G. A., and Johnson, D. E. (2009). *Analysis of Messy Data,
Volume 1: Designed Experiments* (2nd ed.). Chapman and Hall/CRC.
<https://doi.org/10.1201/EBK1584883340>

``` r

fit <- lmm(
  multicentre,
  Y ~ Drug * Time,
  random = ~ 1 | Center,
  repeated = ~ Time | Center:Drug:Subject,
  structure = "ar1",
  method = "REML",
  ddf = "satterthwaite"
)
```

The stored PROC MIXED targets correspond to this model specification.
The validation vignette includes the executable R comparisons, their
observed absolute differences, the exact PROC MIXED program, and every
stored numerical target. The repository does not contain a SAS execution
transcript, so these targets are described as regression benchmarks
rather than a new independent SAS run. Estimates and standard errors use
absolute tolerance `1e-3`; denominator degrees of freedom are compared
at their stored decimal precision.

## Theoretical foundations

The implementation is grounded in the following sources:

- Patterson, H. D., and Thompson, R. (1971), for restricted maximum
  likelihood. <https://doi.org/10.1093/biomet/58.3.545>
- Harville, D. A. (1977), for ML and REML variance-component estimation
  and mixed-model prediction.
  <https://doi.org/10.1080/01621459.1977.10480998>
- LaMotte, L. R. (2007), for a direct derivation of the REML likelihood.
  <https://doi.org/10.1007/s00362-006-0335-6>
- Satterthwaite, F. E. (1946), for approximate denominator degrees of
  freedom. <https://doi.org/10.2307/3002019>
- Fai, A. H.-T., and Cornelius, P. L. (1996), for
  multi-degree-of-freedom approximate F-tests.
  <https://doi.org/10.1080/00949659608811740>
- Searle, S. R., Speed, F. M., and Milliken, G. A. (1980), for
  population marginal means.
  <https://doi.org/10.1080/00031305.1980.10483031>
- Pinheiro, J. C., and Bates, D. M. (2000), for structured covariance
  models. <https://doi.org/10.1007/b98882>

## Boundaries of the evidence

The validation does not cover Kenward-Roger inference, multiple
random-effect terms, generalized responses, likelihood-ratio model
comparison, or large-scale sparse performance because those capabilities
are not implemented in version `0.1.0`.
