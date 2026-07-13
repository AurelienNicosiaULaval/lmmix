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
vcov(object, ...)

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
sigma(object, ...)

# S3 method for class 'lmm'
predict(
  object,
  newdata = NULL,
  re.form = NULL,
  na.action = stats::na.pass,
  ...
)

# S3 method for class 'lmm'
model.matrix(object, data = NULL, ...)

# S3 method for class 'lmm'
model.frame(formula, ...)

# S3 method for class 'lmm'
formula(x, fixed.only = TRUE, ...)

# S3 method for class 'lmm'
terms(x, ...)

# S3 method for class 'lmm'
anova(object, ..., type = 3)

# S3 method for class 'lmm'
confint(object, parm = names(object$coefficients), level = 0.95, ...)
```

## Arguments

- x, object, formula:

  An `lmm` object. The `formula` name follows the argument name of the
  [`model.frame()`](https://rdrr.io/r/stats/model.frame.html) generic.

- ...:

  Additional arguments. Supplying another model to
  [`anova()`](https://rdrr.io/r/stats/anova.html) is not supported.

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

- newdata:

  Optional data frame used for prediction.

- re.form:

  `NULL` includes empirical random effects. Any non-`NULL` value
  requests a fixed-effects-only prediction.

- na.action:

  Missing-value action retained for prediction-method compatibility.

- data:

  Optional data frame used to construct a fixed-effects model matrix.

- fixed.only:

  Whether [`formula()`](https://rdrr.io/r/stats/formula.html) returns
  only the fixed formula. If `FALSE`, it returns fixed, random, and
  repeated formulas in a list.

- parm:

  Fixed-effect parameters requested from
  [`confint()`](https://rdrr.io/r/stats/confint.html), supplied by name
  or position.

- level:

  Confidence level for fixed-effect intervals.

## Value

The return value follows the corresponding base R generic.

## Details

Conditional fitted values and predictions include empirical random
effects. Marginal values use fixed effects only. In prediction from new
data, known grouping levels use their fitted random effect and new
levels receive a random contribution of zero.
