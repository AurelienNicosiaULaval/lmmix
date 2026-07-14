# Print an `lmmix` result table

Printed column headings use human-readable labels without dots. The
underlying object retains its programmatic column names for
compatibility with `broom` and downstream R code.

## Usage

``` r
# S3 method for class 'lmm_table'
print(x, ...)
```

## Arguments

- x:

  An `lmmix` result table.

- ...:

  Additional arguments passed to the tibble print method.

## Value

`x`, invisibly.
