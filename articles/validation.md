# Validation strategy and numerical benchmarks

## Validation principles

Validation is divided into three layers:

1.  structural tests verify positive definiteness, dimensions, parameter
    transformations and input validation;
2.  overlapping model classes are compared with independent R
    implementations;
3.  the combined random-effect and correlated-residual model is checked
    against stored numerical targets for a fully specified PROC MIXED
    analysis.

This separation matters because no single external R implementation
covers every model and inference component supported by `lmmix`.

## Cross-implementation comparisons

The executable test suite includes the following comparisons.

| `lmmix` model | Reference implementation | Quantities compared | Main tolerance |
|:---|:---|:---|:---|
| Random-intercept REML | [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) | Fixed effects, covariance parameters, log-likelihood | `1e-5` or smaller |
| Random-slope REML | [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) | Fixed effects, random covariance, log-likelihood | `1e-4` or smaller |
| Random-intercept ML | [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) | Fixed effects and log-likelihood | `1e-6` |
| Marginal CS and AR(1) | [`nlme::gls()`](https://rdrr.io/pkg/nlme/man/gls.html) | Fixed effects, residual parameters, log-likelihood | `1e-5` |
| Marginal AR(1) | [`mmrm::mmrm()`](https://openpharma.github.io/mmrm/latest-tag/reference/mmrm.html) | Fixed effects, full residual covariance, log-likelihood | `1e-5` |
| Random-intercept Satterthwaite inference | [`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html) | Estimates, standard errors and coefficient degrees of freedom | `1e-3` or smaller |
| Type III Satterthwaite tests | [`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html) | Numerator df, denominator df and F statistics | `1e-3` or smaller |

These comparisons are restricted to cases where the covariance and
inference definitions overlap. Agreement in one model class is not
presented as evidence for an unsupported feature.

## Combined covariance structures

The included multilocation data are used to fit random center effects
together with each supported residual structure: `id`, `cs`, `ar1`,
`toep`, and `un`. For every fit, tests require optimizer convergence, a
positive-definite likelihood Hessian and a positive-definite marginal
covariance matrix.

The data reproduce the repeated-measures example in Section 28.3 of
Milliken and Johnson ([2009](#ref-milliken2009)). They contain 153 rows,
including 28 missing responses. The complete-case model therefore uses
125 observations.

``` r

library(lmmix)

fit <- lmm(
  data = multicentre,
  formula = Y ~ Drug * Time,
  random = ~ 1 | Center,
  repeated = ~ Time | Center:Drug:Subject,
  structure = "ar1",
  method = "REML",
  ddf = "satterthwaite"
)

fit$convergence[c("code", "hessian_positive_definite")]
#> $code
#> [1] 0
#> 
#> $hessian_positive_definite
#> [1] TRUE
```

## Stored PROC MIXED regression targets

The combined AR(1) regression targets correspond to the following
program. The RANDOM and REPEATED statements represent the two covariance
sources, and `DDFM=SATTERTHWAITE` requests row-specific denominator
degrees of freedom as documented by SAS Institute Inc.
([2015](#ref-sas2015)).

``` sas
proc mixed data=multicentre method=reml;
  class Center Drug Subject Time;
  model Y = Drug Time Drug*Time / ddfm=satterth solution;
  random Center;
  repeated Time / type=ar(1) subject=Subject(Center*Drug);
  lsmeans Drug Time / diff;
run;
```

The repository stores numerical targets for this specification. It does
not include a SAS execution transcript. The program above is provided so
the benchmark can be rerun in a licensed SAS environment.

| Quantity           | Level or contrast | Estimate | Standard error |   df |
|:-------------------|:------------------|---------:|---------------:|-----:|
| Drug marginal mean | 1                 |  13.1568 |         1.5290 | 2.98 |
| Drug marginal mean | 2                 |  18.1484 |         1.5324 | 3.01 |
| Drug marginal mean | 3                 |  17.0516 |         1.5327 | 3.01 |
| Time marginal mean | 1                 |  14.9543 |         1.3961 | 2.08 |
| Time marginal mean | 2                 |  15.6700 |         1.3985 | 2.09 |
| Time marginal mean | 3                 |  17.7324 |         1.4032 | 2.12 |
| Drug difference    | 1 minus 2         |  -4.9916 |         1.0990 | 46.1 |

The remaining stored pairwise targets are:

| Factor | Contrast  | Estimate | Standard error |   df |      p-value |
|:-------|:----------|---------:|---------------:|-----:|-------------:|
| Drug   | 1 minus 3 |  -3.8948 |         1.0987 | 46.1 |       0.0009 |
| Drug   | 2 minus 3 |   1.0968 |         1.1043 | 46.9 |       0.3257 |
| Time   | 1 minus 2 |  -0.7157 |         0.1845 | 68.1 |       0.0002 |
| Time   | 1 minus 3 |  -2.7781 |         0.2714 | 73.2 | below 0.0001 |
| Time   | 2 minus 3 |  -2.0624 |         0.2056 | 68.6 | below 0.0001 |

Estimates and standard errors are compared with absolute tolerance
`1e-3`. Degrees of freedom are compared at the decimal precision stored
with each target.

## What the validation does not establish

The test suite does not validate Kenward-Roger inference because that
method is not implemented. It does not establish large-sample
scalability because the current likelihood uses dense covariance
matrices. The stored PROC MIXED targets are regression tests, not a
substitute for a new independent SAS run. These boundaries are stated
explicitly so the validation claims do not exceed the available
evidence.

Milliken, George A., and Dallas E. Johnson. 2009. *Analysis of Messy
Data, Volume 1: Designed Experiments*. 2nd ed. Chapman; Hall/CRC.
<https://doi.org/10.1201/EBK1584883340>.

SAS Institute Inc. 2015. *SAS/STAT 14.1 User’s Guide: The MIXED
Procedure*. Cary, NC.
<https://support.sas.com/documentation/cdl/en/statug/68162/HTML/default/statug_mixed_syntax.htm>.
