# Support for `emmeans`

These methods allow
[`emmeans::emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.html)
to recover the fixed-effects reference grid and use model-specific
Satterthwaite degrees of freedom.

## Usage

``` r
# S3 method for class 'lmm'
recover_data(object, frame = object$model_frame, ...)

# S3 method for class 'lmm'
emm_basis(object, trms, xlev, grid, ...)
```

## Arguments

- object:

  An `lmm` object.

- frame:

  Model frame.

- ...:

  Additional arguments passed by `emmeans`.

- trms, xlev, grid:

  Arguments supplied by `emmeans`.

## Value

An `emmeans` data-recovery object.

A basis list used by `emmeans`.
