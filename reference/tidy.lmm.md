# Tidy an `lmm` model

Tidy an `lmm` model

## Usage

``` r
# S3 method for class 'lmm'
tidy(
  x,
  effects = c("fixed", "ran_pars", "all"),
  conf.int = FALSE,
  conf.level = 0.95,
  ...
)
```

## Arguments

- x:

  An `lmm` object.

- effects:

  Effects to return: fixed effects, covariance parameters, or both.

- conf.int:

  Whether to include confidence limits for fixed effects.

- conf.level:

  Confidence level.

- ...:

  Additional arguments.

## Value

A tibble.
