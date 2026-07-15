# Standard methods for `lmm` fits

Methods are provided for extracting likelihood quantities, fitted
values, residuals, model matrices and formulas, as well as for
prediction, type III testing and fixed-effect confidence intervals.

## Usage

``` r
# S3 method for class 'lmm'
print(x, ...)

# S3 method for class 'lmm'
coef(object, ...)

# S3 method for class 'lmm'
vcov(object, adjusted = TRUE, ...)

# S3 method for class 'lmm'
logLik(object, ...)

# S3 method for class 'lmm'
AIC(object, ..., k = 2)

# S3 method for class 'lmm'
BIC(object, ...)

# S3 method for class 'lmm'
deviance(object, ...)

# S3 method for class 'lmm'
nobs(object, ...)

# S3 method for class 'lmm'
fitted(object, type = c("conditional", "marginal"), ...)

# S3 method for class 'lmm'
residuals(object, type = c("response", "pearson", "marginal"), ...)

# S3 method for class 'lmm'
plot(x, which = c("residuals", "qq", "fitted"), ...)

# S3 method for class 'lmm'
sigma(object, ...)

# S3 method for class 'lmm'
predict(
  object,
  newdata = NULL,
  re.form = NULL,
  allow.new.levels = FALSE,
  se.fit = FALSE,
  interval = c("none", "confidence", "prediction"),
  level = 0.95,
  na.action = stats::na.pass,
  ...
)

# S3 method for class 'lmm'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'lmm'
update(object, formula., ..., evaluate = TRUE)

# S3 method for class 'lmm'
model.matrix(object, data = NULL, ...)

# S3 method for class 'lmm'
model.frame(formula, ...)

# S3 method for class 'lmm'
formula(x, fixed.only = TRUE, ...)

# S3 method for class 'lmm'
terms(x, ...)

# S3 method for class 'lmm'
anova(
  object,
  ...,
  type = 3,
  refit = TRUE,
  test = c("chisq", "parametric.bootstrap"),
  nsim = 199L,
  seed = NULL
)

# S3 method for class 'lmm'
confint(object, parm = names(object$coefficients), level = 0.95, ...)
```

## Arguments

- x, object, formula:

  An `lmm` object. The `formula` name follows the argument name of the
  [`model.frame()`](https://rdrr.io/r/stats/model.frame.html) generic.

- ...:

  Additional fitted `lmm` models for likelihood-ratio comparison, or
  arguments passed to the corresponding method.

- adjusted:

  Whether [`vcov()`](https://rdrr.io/r/stats/vcov.html) returns the
  Kenward-Roger-adjusted covariance when available.

- k:

  Penalty per parameter used by
  [`AIC()`](https://rdrr.io/r/stats/AIC.html).

- type:

  For [`fitted()`](https://rdrr.io/r/stats/fitted.values.html), either
  `"conditional"` or `"marginal"`. For
  [`residuals()`](https://rdrr.io/r/stats/residuals.html), one of
  `"response"`, `"pearson"`, or `"marginal"`. For
  [`anova()`](https://rdrr.io/r/stats/anova.html), only type `3` is
  implemented.

- which:

  Diagnostic plot to draw: standardized residuals against fitted values,
  a normal quantile-quantile plot, or observed against fitted values.

- newdata:

  Optional data frame used for prediction.

- re.form:

  `NULL` includes every empirical random-effect term, `NA` or `~ 0`
  excludes all random effects, and a random-effects formula or list of
  formulas selects a subset of fitted terms.

- allow.new.levels:

  Whether conditional predictions may include new grouping levels. Their
  random contribution is zero when allowed.

- se.fit:

  Whether [`predict()`](https://rdrr.io/r/stats/predict.html) returns
  standard errors for the expected response. These standard errors use
  the fixed-effect covariance and condition on any empirical random
  effects included through `re.form`.

- interval:

  Prediction interval type. `"confidence"` describes the expected
  response. `"prediction"` additionally includes the fitted residual
  variance but not uncertainty in empirical random effects.

- level:

  Confidence level for fixed-effect and covariance intervals.

- na.action:

  Missing-value action retained for prediction-method compatibility.

- nsim:

  Number of simulations for a parametric-bootstrap comparison.

- seed:

  Optional integer seed for a reproducible parametric bootstrap.

- formula.:

  Optional formula update accepted by
  [`update()`](https://rdrr.io/r/stats/update.html).

- evaluate:

  Whether [`update()`](https://rdrr.io/r/stats/update.html) evaluates
  the updated call.

- data:

  Optional data frame used to construct a fixed-effects model matrix.

- fixed.only:

  Whether [`formula()`](https://rdrr.io/r/stats/formula.html) returns
  only the fixed formula. If `FALSE`, it returns fixed, random, and
  repeated formulas in a list.

- refit:

  Whether REML models that differ in fixed effects are refitted with ML
  before comparison.

- test:

  Reference test for likelihood-ratio model comparisons. The default
  `"chisq"` uses the usual asymptotic chi-squared reference. The
  `"parametric.bootstrap"` option simulates under each smaller model and
  is appropriate when covariance parameters can lie on a boundary.

- parm:

  Parameters requested from
  [`confint()`](https://rdrr.io/r/stats/confint.html), supplied by name
  or fixed effect position. Use `"beta_"` for every fixed effect and
  `"theta_"` for every covariance parameter.

## Value

The return value follows the corresponding base R generic.

## Details

Conditional fitted values and predictions include empirical random
effects. Marginal values use fixed effects only. In prediction from new
data, known grouping levels use their fitted random effect. New levels
are rejected by default and receive a random contribution of zero only
when explicitly allowed.

[`ranef()`](https://aureliennicosiaulaval.github.io/lmmix/reference/ranef.md)
always returns a named list with one tibble per random-effect term,
including for models with a single term. It returns an empty list for
marginal models. [`sigma()`](https://rdrr.io/r/stats/sigma.html) is
defined as `sqrt(mean(diag(R)))`; it is a descriptive residual scale
when residual variances differ across observations.

## References

Self, S. G., and Liang, K.-Y. (1987). Asymptotic properties of maximum
likelihood estimators and likelihood ratio tests under nonstandard
conditions. *Journal of the American Statistical Association*, 82(398),
605-610.
[doi:10.1080/01621459.1987.10478472](https://doi.org/10.1080/01621459.1987.10478472)

Davison, A. C., and Hinkley, D. V. (1997). *Bootstrap Methods and Their
Application*. Cambridge University Press.
[doi:10.1017/CBO9780511802843](https://doi.org/10.1017/CBO9780511802843)
