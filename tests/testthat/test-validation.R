test_that("input validation gives informative errors", {
  expect_error(lmm(list(), y ~ x), "data frame or tibble")
  expect_error(lmm(data.frame(x = 1:3), y ~ x), "Unknown variable")
  expect_error(
    lmm(data.frame(y = 1:3, x = 1:3), y ~ x, structure = "ar1"),
    "repeated"
  )
  expect_error(
    lmm(
      data.frame(y = 1:4, time = rep(1:2, 2), id = 1),
      y ~ 1,
      repeated = ~ time | id,
      structure = "ar1"
    ),
    "at most one observation"
  )
})

test_that("data can be piped into lmm", {
  fit <- orthodont_data() |>
    lmm(formula = distance ~ age + Sex, random = ~ 1 | Subject)
  expect_s3_class(fit, "lmm")
})

test_that("the optim backend is available", {
  fit <- lmm(
    orthodont_data(),
    distance ~ age + Sex,
    random = ~ 1 | Subject,
    control = lmm_control(
      optimizer = "optim",
      optim_method = "BFGS",
      max_iter = 500
    )
  )

  expect_identical(fit$convergence$optimizer, "optim")
  expect_identical(fit$convergence$code, 0L)
})

test_that("optimization controls are validated", {
  expect_error(lmm_control(max_iter = 0), "positive")
  expect_error(lmm_control(rel_tol = 0), "positive")
  expect_error(lmm_control(x_tol = 0), "positive")
  expect_error(lmm_control(lower = 1, upper = 0), "smaller")
  expect_error(lmm_control(initial = "bad"), "numeric vector")
  expect_error(
    lmm(
      orthodont_data(),
      distance ~ age + Sex,
      random = ~ 1 | Subject,
      control = lmm_control(initial = 1:3)
    ),
    "finite numeric"
  )
})
