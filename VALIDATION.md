# Validation scope

## Evidence layers

The package uses four complementary validation layers:

1. structural tests for covariance parameterization, positive definiteness,
   dimensions, missing-data handling, input validation, methods, and output
   classes;
2. numerical comparisons with independent R implementations in overlapping
   model classes;
3. published targets from compatible official PROC MIXED examples;
4. stored numerical regression targets for a combined random-effect and AR(1)
   residual model specified in PROC MIXED.

No single external R package covers every covariance and inference component
implemented by `lmmix`. Claims are therefore restricted to the quantities
tested in each overlapping case.

## Independent R comparisons

| Model | Reference | Quantities checked | Tolerance |
|:--|:--|:--|:--|
| Random-intercept REML | `nlme::lme()` | Fixed effects and log-likelihood | absolute `1e-6` |
| Random-intercept REML | `nlme::lme()` | Covariance parameters | absolute `1e-4` |
| Random-slope REML | `nlme::lme()` | Fixed effects and log-likelihood | absolute `1e-5` |
| Random-slope REML | `nlme::lme()` | Random covariance and residual variance | absolute `1e-4` |
| Random-intercept ML | `nlme::lme()` | Fixed effects and log-likelihood | `1e-6` |
| Marginal CS and AR(1) | `nlme::gls()` | Fixed effects, residual parameters, log-likelihood | `1e-5` |
| Marginal AR(1) | `mmrm::mmrm()` | Fixed effects and log-likelihood | absolute `1e-5` |
| Marginal AR(1) | `mmrm::mmrm()` | Residual covariance matrix | absolute `1e-4` |
| Satterthwaite coefficients | `lmerTest::lmer()` | Estimates and standard errors | absolute `1e-5` |
| Satterthwaite coefficients | `lmerTest::lmer()` | Denominator df | absolute `5e-3` |
| Type III Satterthwaite tests | `lmerTest::lmer()` | Numerator df, denominator df, F statistics | absolute `1e-6`, `5e-3`, and `5e-3` |
| Marginal AR(1) Kenward-Roger | `mmrm::mmrm()` | Fixed effects, adjusted standard errors, denominator df | absolute `1e-5`, `3e-3`, and `5e-3` |
| Crossed random intercepts | `lme4::lmer()` | Log-likelihood and covariance parameters | absolute `1e-5` and `1e-4` |
| Nested fixed-effect models | `lme4::anova()` | Likelihood-ratio statistic and p-value | absolute `1e-5` |

The combined-model tests fit `id`, `cs`, `ar1`, `toep`, and `un` residual
structures together with a random center intercept. Every fit must converge
with a positive-definite likelihood Hessian and marginal covariance matrix.

## Official PROC MIXED examples

Four examples from the official SAS/STAT 14.3 documentation overlap exactly
with the current `lmmix` model space.

| SAS example | lmmix model | Quantities checked |
|:--|:--|:--|
| 79.1 Split-Plot Design | Two independent random-intercept terms, REML | Variance components, restricted likelihood criterion, type III tests |
| 79.2 Repeated Measures | ML with UN and CS residual covariance | Covariance matrices, likelihood criterion, fixed-effect solutions, standard errors, type III statistics |
| 79.5 Random Coefficients | Correlated random intercept and slope, REML | Random covariance, residual variance, likelihood criterion, fixed effects, BLUPs, type III test |
| 79.6 Line-Source Sprinkler | Three random-intercept terms and a four-band Toeplitz residual covariance, REML | Random variances, residual covariance bands, restricted likelihood criterion |

The machine-readable targets and the runnable SAS program are installed under
`validation/sas`. They come from the published output tables in SAS Institute
Inc. (2017), *SAS/STAT 14.3 User's Guide: The MIXED Procedure*:
<https://support.sas.com/documentation/onlinedoc/stat/>.

These are documentation-based comparisons. No claim is made that SAS was run
inside this repository. The full SAS manual is not redistributed.

## Multilocation data

The included `multicentre` data reproduce the multilocation repeated-measures
example in Section 28.3 of Milliken and Johnson (2009). The data contain 153
rows, including 28 missing responses, so the model below uses 125 complete
observations.

Milliken, G. A., and Johnson, D. E. (2009). *Analysis of Messy Data, Volume 1:
Designed Experiments* (2nd ed.). Chapman and Hall/CRC.
<https://doi.org/10.1201/EBK1584883340>

```r
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

The stored PROC MIXED targets correspond to this model specification. The
validation vignette includes the executable R comparisons, their observed
absolute differences, the exact PROC MIXED program, and every stored numerical
target. The repository does not contain a SAS execution transcript, so these
targets are described as regression benchmarks rather than a new independent
SAS run. Estimates and standard errors use absolute tolerance `1e-3`;
denominator degrees of freedom are compared at their stored decimal precision.

## Theoretical foundations

The implementation is grounded in the following sources:

* Patterson, H. D., and Thompson, R. (1971), for restricted maximum
  likelihood. <https://doi.org/10.1093/biomet/58.3.545>
* Harville, D. A. (1977), for ML and REML variance-component estimation and
  mixed-model prediction. <https://doi.org/10.1080/01621459.1977.10480998>
* LaMotte, L. R. (2007), for a direct derivation of the REML likelihood.
  <https://doi.org/10.1007/s00362-006-0335-6>
* Satterthwaite, F. E. (1946), for approximate denominator degrees of freedom.
  <https://doi.org/10.2307/3002019>
* Fai, A. H.-T., and Cornelius, P. L. (1996), for multi-degree-of-freedom
  approximate F-tests. <https://doi.org/10.1080/00949659608811740>
* Kenward, M. G., and Roger, J. H. (1997), for small-sample covariance and
  denominator-df adjustment. <https://doi.org/10.2307/2533558>
* Self, S. G., and Liang, K.-Y. (1987), for likelihood-ratio tests with
  parameters on the boundary. <https://doi.org/10.1080/01621459.1987.10478472>
* Davison, A. C., and Hinkley, D. V. (1997), for parametric-bootstrap
  inference. <https://doi.org/10.1017/CBO9780511802843>
* Searle, S. R., Speed, F. M., and Milliken, G. A. (1980), for population
  marginal means. <https://doi.org/10.1080/00031305.1980.10483031>
* Pinheiro, J. C., and Bates, D. M. (2000), for structured covariance models.
  <https://doi.org/10.1007/b98882>

## Boundaries of the evidence

The validation covers Kenward-Roger inference in overlapping random-effects
and marginal models, multiple crossed random intercepts, and fixed-effect
likelihood-ratio comparisons. Generalized responses and large-scale sparse
performance remain outside the package scope. The line-source example supplies
an official target for a combined random-effect and correlated-residual model.
Combined-model Kenward-Roger inference has structural tests but no fresh
independent SAS execution.
