# Control numerical optimization for `lmm()`

Control numerical optimization for [`lmm()`](lmm.md)

## Usage

``` r
lmm_control(
  optimizer = c("nlminb", "optim"),
  optim_method = "BFGS",
  max_iter = 1000L,
  rel_tol = 1e-08,
  x_tol = 1e-08,
  initial = NULL,
  lower = -20,
  upper = 20,
  deriv_method = "Richardson"
)
```

## Arguments

- optimizer:

  Optimizer to use, either `"nlminb"` or `"optim"`.

- optim_method:

  Method passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) when
  `optimizer = "optim"`.

- max_iter:

  Maximum number of optimizer iterations.

- rel_tol:

  Relative convergence tolerance.

- x_tol:

  Parameter convergence tolerance for
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- initial:

  Optional numeric vector of unconstrained starting values.

- lower, upper:

  Bounds for unconstrained covariance parameters.

- deriv_method:

  Numerical differentiation method passed to `numDeriv`.

## Value

An object of class `lmm_control`.
