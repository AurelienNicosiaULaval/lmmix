# lmmix

`lmmix` fits Gaussian linear mixed models by explicitly optimizing the ML or
REML criterion for

\[
V = ZGZ^\mathsf{T} + R.
\]

The package is designed for models that combine random effects with correlated
residuals. It provides Satterthwaite denominator degrees of freedom, type III
tests, estimated marginal means, pairwise contrasts, broom methods, and
`emmeans` interoperability.

The package name is `lmmix`; the main fitting function is `lmm()`.

## Installation

Install the package from a local checkout:

```r
devtools::install("path/to/lmmix")
```

## Combined random-effect and residual-correlation model

The included `multicentre` data are simulated. They reproduce the design in
the development brief, not the original observations used for the SAS
reference table.

```r
library(lmmix)

fit <- lmm(
  data = multicentre,
  formula = Y ~ Drug * Time,
  random = ~ 1 | Center,
  repeated = ~ Time | Center:Subject,
  structure = "ar1",
  method = "REML",
  ddf = "satterthwaite"
)

summary(fit)
anova(fit)
lsmeans(fit, pairwise ~ Drug)
```

The data-first interface also works with the base R pipe:

```r
fit <- multicentre |>
  lmm(
    formula = Y ~ Drug * Time,
    random = ~ 1 | Center,
    repeated = ~ Time | Center:Subject,
    structure = "ar1"
  )
```

## Supported covariance structures

The residual structures are:

* `id`: homogeneous independent residuals;
* `cs`: compound symmetry;
* `ar1`: first-order autoregression, including negative correlation;
* `toep`: positive-definite Toeplitz covariance parameterized through partial
  autocorrelations;
* `un`: unstructured covariance parameterized through a Cholesky factor.

Random intercepts and random slopes use an unstructured Cholesky-parameterized
`G` matrix.

## Validation status

Executable tests compare the package with `nlme` for ML, REML, random
intercepts, random slopes, compound symmetry, and AR(1). Satterthwaite tests
are compared with `lmerTest` where the models overlap.

The original multi-site observations needed to verify the numerical SAS
targets were not supplied. Those target tests are included and explicitly
skipped. They must not be treated as passed until the original data are added.

## References

* Satterthwaite, F. E. (1946). An approximate distribution of estimates of
  variance components. *Biometrics Bulletin*, 2(6), 110-114.
  <https://doi.org/10.2307/3002019>
* Pinheiro, J. C., and Bates, D. M. (2000). *Mixed-Effects Models in S and
  S-PLUS*. Springer. <https://doi.org/10.1007/b98882>
