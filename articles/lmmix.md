# Getting started with lmmix

`lmmix` fits the Gaussian model

``` math
Y = X\beta + Zu + \varepsilon,
\qquad
V = ZGZ^\mathsf{T} + R,
```

by directly optimizing a profiled ML or REML objective. The model
formulation and covariance decomposition follow standard mixed-model
theory ([Harville 1977](#ref-harville1977); [Pinheiro and Bates
2000](#ref-pinheiro2000)).

``` r

library(lmmix)
```

## A model with two covariance sources

The example data contain a random center effect and repeated
observations within subjects.

``` r

fit <- lmm(
  multicentre,
  Y ~ Drug * Time,
  random = ~ 1 | Center,
  repeated = ~ Time | Center:Drug:Subject,
  structure = "ar1"
)

fit
#> Linear mixed model fit by REML
#> Formula: `Y ~ Drug * Time`
#> Random: `~1 | Center`
#> Repeated: `~Time | Center:Drug:Subject` (AR1)
#> Log-likelihood: -245.313
#> Convergence code: 0
```

The summary contains coefficient-level tests, type III tests, covariance
parameters, and information criteria.

``` r

summary(fit)
#> 
#> ── Linear mixed model ──────────────────────────────────────────────────────────
#> Estimation: REML
#> Denominator df: satterthwaite
#> Residual covariance: AR1
#> 
#> ── Fixed effects ──
#> 
#> # A tibble: 9 × 8
#>   term        estimate std.error statistic    df     p.value conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl> <dbl>       <dbl>    <dbl>     <dbl>
#> 1 (Intercept)   12.1       1.54      7.83   3.06 0.00405        7.21     16.9  
#> 2 Drug2          4.76      1.12      4.25  49.8  0.0000928      2.51      7.02 
#> 3 Drug3          3.94      1.12      3.52  49.8  0.000939       1.69      6.19 
#> 4 Time2          0.921     0.294     3.13  67.5  0.00255        0.335     1.51 
#> 5 Time3          2.39      0.441     5.42  72.1  0.000000750    1.51      3.27 
#> 6 Drug2:Time2   -0.143     0.438    -0.326 67.9  0.746         -1.02      0.732
#> 7 Drug3:Time2   -0.474     0.448    -1.06  68.0  0.294         -1.37      0.421
#> 8 Drug2:Time3    0.823     0.666     1.24  72.7  0.221         -0.505     2.15 
#> 9 Drug3:Time3    0.334     0.644     0.519 73.1  0.605         -0.949     1.62
#> ── Type III tests ──
#> # A tibble: 3 × 5
#>   term      num.df den.df statistic  p.value
#>   <chr>      <dbl>  <dbl>     <dbl>    <dbl>
#> 1 Drug           2   46.4     11.4  9.23e- 5
#> 2 Time           2   68.5     59.3  1.16e-15
#> 3 Drug:Time      4   68.4      1.35 2.60e- 1
#> ── Covariance parameters ──
#> # A tibble: 3 × 5
#>   group    term        component estimate std.error
#>   <chr>    <chr>       <chr>        <dbl>     <dbl>
#> 1 Center   (Intercept) var          5.17     5.73  
#> 2 Residual ar1         var         10.7      2.14  
#> 3 Residual ar1         cor          0.935    0.0169
#> ── Information criteria ──
#>    logLik       AIC       BIC  deviance 
#> -245.3129  514.6257  548.5655  490.6257
```

## Estimated marginal means

Nuisance-factor levels receive equal weights. Numeric covariates not
requested in `specs` are held at their observed means. This construction
targets population marginal means in the sense of Searle et al.
([1980](#ref-searle1980)).

``` r

lsmeans(fit, ~Drug)
#> # A tibble: 3 × 8
#>   Drug  estimate std.error    df statistic p.value conf.low conf.high
#>   <fct>    <dbl>     <dbl> <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#> 1 1         13.2      1.53  2.98      8.60 0.00339     8.27      18.0
#> 2 2         18.1      1.53  3.01     11.8  0.00128    13.3       23.0
#> 3 3         17.1      1.53  3.01     11.1  0.00154    12.2       21.9
lsmeans(fit, pairwise ~ Drug)
#> ── Estimated marginal means ──
#> 
#> # A tibble: 3 × 8
#>   Drug  estimate std.error    df statistic p.value conf.low conf.high
#>   <fct>    <dbl>     <dbl> <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#> 1 1         13.2      1.53  2.98      8.60 0.00339     8.27      18.0
#> 2 2         18.1      1.53  3.01     11.8  0.00128    13.3       23.0
#> 3 3         17.1      1.53  3.01     11.1  0.00154    12.2       21.9
#> ── Pairwise contrasts ──
#> # A tibble: 3 × 8
#>   contrast estimate std.error    df statistic   p.value conf.low conf.high
#>   <chr>       <dbl>     <dbl> <dbl>     <dbl>     <dbl>    <dbl>     <dbl>
#> 1 1 - 2       -4.99      1.10  46.1    -4.54  0.0000400    -7.20     -2.78
#> 2 1 - 3       -3.89      1.10  46.1    -3.54  0.000914     -6.11     -1.68
#> 3 2 - 3        1.10      1.10  46.9     0.993 0.326        -1.12      3.32
```

## Tidy output

``` r

generics::tidy(fit)
#> # A tibble: 9 × 7
#>   effect term        estimate std.error statistic    df     p.value
#>   <chr>  <chr>          <dbl>     <dbl>     <dbl> <dbl>       <dbl>
#> 1 fixed  (Intercept)   12.1       1.54      7.83   3.06 0.00405    
#> 2 fixed  Drug2          4.76      1.12      4.25  49.8  0.0000928  
#> 3 fixed  Drug3          3.94      1.12      3.52  49.8  0.000939   
#> 4 fixed  Time2          0.921     0.294     3.13  67.5  0.00255    
#> 5 fixed  Time3          2.39      0.441     5.42  72.1  0.000000750
#> 6 fixed  Drug2:Time2   -0.143     0.438    -0.326 67.9  0.746      
#> 7 fixed  Drug3:Time2   -0.474     0.448    -1.06  68.0  0.294      
#> 8 fixed  Drug2:Time3    0.823     0.666     1.24  72.7  0.221      
#> 9 fixed  Drug3:Time3    0.334     0.644     0.519 73.1  0.605
generics::glance(fit)
#> # A tibble: 1 × 8
#>   logLik   AIC   BIC deviance    df  nobs convergence method
#>    <dbl> <dbl> <dbl>    <dbl> <int> <int>       <int> <chr> 
#> 1  -245.  515.  549.     491.    12   125           0 REML
head(generics::augment(fit))
#> # A tibble: 6 × 8
#>   Center Drug  Subject Time      Y .fitted .resid .std.resid
#>   <fct>  <fct> <fct>   <fct> <dbl>   <dbl>  <dbl>      <dbl>
#> 1 R      1     1       1        17    13.0  4.05      1.24  
#> 2 R      1     2       1        12    13.0 -0.954    -0.292 
#> 3 R      1     2       2        14    13.9  0.125     0.0383
#> 4 R      1     2       3        15    15.3 -0.346    -0.106 
#> 5 R      1     3       1        12    13.0 -0.954    -0.292 
#> 6 R      1     3       2        11    13.9 -2.87     -0.880
```

When `emmeans` is installed, the package supplies the required basis and
uses the model’s Satterthwaite degrees of freedom.

``` r

emmeans::emmeans(fit, ~Drug)
#> NOTE: Results may be misleading due to involvement in interactions
#>  Drug emmean   SE   df lower.CL upper.CL
#>  1      13.2 1.53 2.98     8.27     18.0
#>  2      18.1 1.53 3.01    13.28     23.0
#>  3      17.1 1.53 3.01    12.18     21.9
#> 
#> Results are averaged over the levels of: Time 
#> Confidence level used: 0.95
```

Harville, David A. 1977. “Maximum Likelihood Approaches to Variance
Component Estimation and to Related Problems.” *Journal of the American
Statistical Association* 72 (358): 320–38.
<https://doi.org/10.1080/01621459.1977.10480998>.

Pinheiro, José C., and Douglas M. Bates. 2000. *Mixed-Effects Models in
s and s-PLUS*. Springer. <https://doi.org/10.1007/b98882>.

Searle, S. R., F. M. Speed, and G. A. Milliken. 1980. “Population
Marginal Means in the Linear Model: An Alternative to Least Squares
Means.” *The American Statistician* 34 (4): 216–21.
<https://doi.org/10.1080/00031305.1980.10483031>.
