# Fit a Gaussian linear mixed model

`lmm()` explicitly evaluates a profiled ML or REML criterion for
`V = Z G Z' + R`. It can combine random effects with an independent,
compound-symmetric, AR(1), Toeplitz, or unstructured residual
covariance.

## Usage

``` r
lmm(
  data,
  formula,
  random = NULL,
  repeated = NULL,
  structure = c("id", "cs", "ar1", "toep", "un"),
  method = c("REML", "ML"),
  ddf = c("satterthwaite", "residual", "kenward-roger"),
  control = lmm_control()
)
```

## Arguments

- data:

  A data frame or tibble. Placing `data` first makes the function
  pipe-friendly.

- formula:

  A two-sided fixed-effects formula.

- random:

  Optional one-sided random-effects formula such as `~ 1 | subject` or
  `~ 1 + time | subject`.

- repeated:

  Optional one-sided repeated-measures formula such as
  `~ time | subject`.

- structure:

  Residual covariance structure: `"id"`, `"cs"`, `"ar1"`, `"toep"`, or
  `"un"`.

- method:

  Estimation method: restricted maximum likelihood (`"REML"`) or maximum
  likelihood (`"ML"`).

- ddf:

  Denominator degrees-of-freedom method. `"satterthwaite"` and
  `"residual"` are implemented. `"kenward-roger"` is reserved for future
  implementation.

- control:

  An object returned by [`lmm_control()`](lmm_control.md).

## Value

An object of class `lmm`.

## References

Harville, D. A. (1977). Maximum likelihood approaches to variance
component estimation and to related problems. *Journal of the American
Statistical Association*, 72(358), 320-338.
<https://doi.org/10.1080/01621459.1977.10480998>

LaMotte, L. R. (2007). A direct derivation of the REML likelihood
function. *Statistical Papers*, 48(2), 321-327.
<https://doi.org/10.1007/s00362-006-0335-6>

Pinheiro, J. C., and Bates, D. M. (2000). *Mixed-Effects Models in S and
S-PLUS*. Springer. <https://doi.org/10.1007/b98882>

## Examples

``` r
fit <- lmm(
  data = nlme::Orthodont,
  formula = distance ~ age + Sex,
  random = ~ 1 | Subject
)
summary(fit)
#> 
#> ── Linear mixed model ──────────────────────────────────────────────────────────
#> Estimation: REML
#> Denominator df: satterthwaite
#> Residual covariance: ID
#> 
#> ── Fixed effects ──
#> 
#> # A tibble: 3 × 8
#>   term        estimate std.error statistic    df  p.value conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl> <dbl>    <dbl>    <dbl>     <dbl>
#> 1 (Intercept)   17.7      0.834      21.2   99.4 1.04e-38   16.1      19.4  
#> 2 age            0.660    0.0616     10.7   80.0 3.95e-17    0.538     0.783
#> 3 SexFemale     -2.32     0.761      -3.05  25.0 5.38e- 3   -3.89     -0.753
#> ── Type III tests ──
#> 
#> # A tibble: 2 × 5
#>   term  num.df den.df statistic  p.value
#>   <chr>  <dbl>  <dbl>     <dbl>    <dbl>
#> 1 age        1   80.0    115.   3.95e-17
#> 2 Sex        1   25.0      9.29 5.38e- 3
#> ── Covariance parameters ──
#> 
#> # A tibble: 2 × 5
#>   group    term        component estimate std.error
#>   <chr>    <chr>       <chr>        <dbl>     <dbl>
#> 1 Subject  (Intercept) var           3.27     1.07 
#> 2 Residual id          var           2.05     0.324
#> ── Information criteria ──
#> 
#>    logLik       AIC       BIC  deviance 
#> -218.7563  447.5125  460.9232  437.5125 
```
