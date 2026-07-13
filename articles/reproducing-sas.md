# Multi-site PROC MIXED validation workflow

## Scope of the available data

The `multicentre` data reproduce the multilocation repeated-measures
example described by Milliken and Johnson ([2009](#ref-milliken2009)) in
Section 28.3. The design has three centers, three drugs, three
measurement times, and 153 rows. There are 28 missing responses, so the
fitted model uses 125 complete responses.

## End-to-end model

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

fit
#> Linear mixed model fit by REML
#> Formula: `Y ~ Drug * Time`
#> Random: `~1 | Center`
#> Repeated: `~Time | Center:Drug:Subject` (AR1)
#> Log-likelihood: -245.313
#> Convergence code: 0
VarCorr(fit)
#> # A tibble: 3 × 5
#>   group    term        component estimate std.error
#>   <chr>    <chr>       <chr>        <dbl>     <dbl>
#> 1 Center   (Intercept) var          5.17     5.73  
#> 2 Residual ar1         var         10.7      2.14  
#> 3 Residual ar1         cor          0.935    0.0169
anova(fit)
#> Type III tests with model degrees of freedom
#> # A tibble: 3 × 5
#>   term      num.df den.df statistic  p.value
#>   <chr>      <dbl>  <dbl>     <dbl>    <dbl>
#> 1 Drug           2   46.4     11.4  9.23e- 5
#> 2 Time           2   68.5     59.3  1.16e-15
#> 3 Drug:Time      4   68.4      1.35 2.60e- 1
```

Treatment and time LS-means and treatment differences are generated
without hand-written contrast rows.

``` r

lsmeans(fit, ~Drug)
#> # A tibble: 3 × 8
#>   Drug  estimate std.error    df statistic p.value conf.low conf.high
#>   <fct>    <dbl>     <dbl> <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#> 1 1         13.2      1.53  2.98      8.60 0.00339     8.27      18.0
#> 2 2         18.1      1.53  3.01     11.8  0.00128    13.3       23.0
#> 3 3         17.1      1.53  3.01     11.1  0.00154    12.2       21.9
lsmeans(fit, ~Time)
#> # A tibble: 3 × 8
#>   Time  estimate std.error    df statistic p.value conf.low conf.high
#>   <fct>    <dbl>     <dbl> <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#> 1 1         15.0      1.40  2.08      10.7 0.00753     9.16      20.8
#> 2 2         15.7      1.40  2.09      11.2 0.00670     9.90      21.4
#> 3 3         17.7      1.40  2.12      12.6 0.00495    12.0       23.5
lsmeans(fit, pairwise ~ Drug)
#> 
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

## PROC MIXED benchmark

The stored reference values correspond to the following PROC MIXED
specification. The RANDOM and REPEATED statements represent the two
covariance sources, and `DDFM=SATTERTHWAITE` requests row-specific
denominator degrees of freedom as documented by SAS Institute Inc.
([2026](#ref-sas2026)).

``` sas
proc mixed data=multicentre method=reml;
  class Center Drug Subject Time;
  model Y = Drug Time Drug*Time / ddfm=satterth solution;
  random Center;
  repeated Time / type=ar(1) subject=Subject(Center*Drug);
  lsmeans Drug Time / diff;
run;
```

The stored SAS reference values, with levels oriented in ascending data
order, are:

| Quantity        | Level or contrast | Estimate | Standard error |   df |
|:----------------|:------------------|---------:|---------------:|-----:|
| Drug LS-mean    | 1                 |  13.1568 |         1.5290 | 2.98 |
| Drug LS-mean    | 2                 |  18.1484 |         1.5324 | 3.01 |
| Drug LS-mean    | 3                 |  17.0516 |         1.5327 | 3.01 |
| Time LS-mean    | 1                 |  14.9543 |         1.3961 | 2.08 |
| Time LS-mean    | 2                 |  15.6700 |         1.3985 | 2.09 |
| Time LS-mean    | 3                 |  17.7324 |         1.4032 | 2.12 |
| Drug difference | 1 minus 2         |  -4.9916 |         1.0990 | 46.1 |

The two-sided p-value for the last contrast is below `0.0001`; the
package computes `3.998e-05` before rounding.

The remaining pairwise references are:

| Factor | Contrast  | Estimate | Standard error |   df |      p-value |
|:-------|:----------|---------:|---------------:|-----:|-------------:|
| Drug   | 1 minus 3 |  -3.8948 |         1.0987 | 46.1 |       0.0009 |
| Drug   | 2 minus 3 |   1.0968 |         1.1043 | 46.9 |       0.3257 |
| Time   | 1 minus 2 |  -0.7157 |         0.1845 | 68.1 |       0.0002 |
| Time   | 1 minus 3 |  -2.7781 |         0.2714 | 73.2 | below 0.0001 |
| Time   | 2 minus 3 |  -2.0624 |         0.2056 | 68.6 | below 0.0001 |

## Numerical precision

The benchmark degrees of freedom are retained at the precision shown by
PROC MIXED. Tests compare estimates and standard errors with absolute
tolerance `1e-3`, and degrees of freedom at the displayed precision. The
implementation is also checked against `nlme`, `lmerTest`, and `mmrm` in
model classes where those packages overlap with `lmmix`.

Milliken, George A., and Dallas E. Johnson. 2009. *Analysis of Messy
Data, Volume 1: Designed Experiments*. 2nd ed. Chapman; Hall/CRC.
<https://doi.org/10.1201/EBK1584883340>.

SAS Institute Inc. 2026. *SAS/STAT User’s Guide: The MIXED Procedure*.
<https://documentation.sas.com/doc/en/statug/latest/statug_mixed_syntax01.htm>.
