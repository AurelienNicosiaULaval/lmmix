# Augment data with fitted values and residuals

Augment data with fitted values and residuals

## Usage

``` r
# S3 method for class 'lmm'
augment(x, data = x$data, newdata = NULL, ...)
```

## Arguments

- x:

  An `lmm` object.

- data:

  Data used for augmentation. Defaults to the analysis data.

- newdata:

  Optional new data. If supplied, it takes precedence over `data`.

- ...:

  Additional arguments passed to
  [`predict()`](https://rdrr.io/r/stats/predict.html).

## Value

A tibble with `.fitted`, `.resid`, and `.std.resid` columns.
