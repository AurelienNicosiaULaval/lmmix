# Augment data with fitted values and residuals

Augment data with fitted values and residuals

## Usage

``` r
# S3 method for class 'lmm'
augment(x, data = NULL, newdata = NULL, ...)
```

## Arguments

- x:

  An `lmm` object.

- data:

  Optional data used for augmentation. By default, excluded rows are
  restored when the model used
  [`stats::na.exclude()`](https://rdrr.io/r/stats/na.fail.html).

- newdata:

  Optional new data. If supplied, it takes precedence over `data`.

- ...:

  Additional arguments passed to
  [`predict()`](https://rdrr.io/r/stats/predict.html).

## Value

A tibble with `.fitted`, `.resid`, and `.std.resid` columns.
