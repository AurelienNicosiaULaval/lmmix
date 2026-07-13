# Validation scope and theoretical sources

## Data provenance

The `multicentre` data reproduce the multilocation repeated-measures example
in Section 28.3 of Milliken and Johnson (2009). The data contain 153 rows,
including 28 missing responses, so the reference model uses 125 observations.

Milliken, G. A., and Johnson, D. E. (2009). *Analysis of Messy Data, Volume 1:
Designed Experiments* (2nd ed.). Chapman and Hall/CRC.
<https://doi.org/10.1201/EBK1584883340>

## Reference model

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

The active tests compare covariance parameters, fixed effects, type III tests,
LS-means, pairwise differences, and rowwise Satterthwaite degrees of freedom
with stored SAS PROC MIXED benchmarks for this specification. The exact SAS
statements are given in the validation vignette. PROC MIXED algorithms and
syntax are documented in the current SAS/STAT User's Guide:
<https://documentation.sas.com/doc/en/statug/latest/statug_mixed_syntax01.htm>.

Estimates and standard errors are tested with absolute tolerance `1e-3`.
Degrees of freedom are compared at the precision stored with the benchmark.
Overlapping model classes are independently compared with `nlme`, `lmerTest`,
and `mmrm`.

## Theoretical foundations

The implementation is grounded in the following primary references:

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
* Searle, S. R., Speed, F. M., and Milliken, G. A. (1980), for population
  marginal means. <https://doi.org/10.1080/00031305.1980.10483031>
* Pinheiro, J. C., and Bates, D. M. (2000), for structured covariance models.
  <https://doi.org/10.1007/b98882>
