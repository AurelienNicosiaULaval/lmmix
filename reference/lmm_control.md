# Control numerical optimization for `lmm()`

Control numerical optimization for
[`lmm()`](https://aureliennicosiaulaval.github.io/lmmix/reference/lmm.md)

## Usage

``` r
lmm_control(
  optimizer = c("auto", "nlminb", "optim"),
  optim_method = "BFGS",
  max_iter = 1000L,
  rel_tol = 1e-08,
  x_tol = 1e-08,
  initial = NULL,
  lower = -20,
  upper = 20,
  deriv_method = "Richardson",
  max_restarts = 3L,
  restart_scale = 0.25
)
```

## Arguments

- optimizer:

  Optimizer strategy: `"auto"`, `"nlminb"`, or `"optim"`.

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

- max_restarts:

  Maximum number of deterministic restart attempts after an unsuccessful
  initial fit.

- restart_scale:

  Size of deterministic perturbations to starting values.

## Value

An object of class `lmm_control`.
