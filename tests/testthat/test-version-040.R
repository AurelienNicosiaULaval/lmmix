test_that("blockwise GLS agrees with the dense calculation", {
  fit <- fit_orthodont_intercept(ddf = "residual")
  expect_gt(length(fit$design$likelihood_blocks), 1L)

  blockwise <- gls_at_eta(fit$eta, fit$design)
  dense_design <- fit$design
  dense_design$likelihood_blocks <- list(seq_len(nobs(fit)))
  dense <- gls_at_eta(fit$eta, dense_design)

  expect_equal(blockwise$beta, dense$beta, tolerance = 1e-10)
  expect_equal(blockwise$beta_vcov, dense$beta_vcov, tolerance = 1e-10)
  expect_equal(blockwise$quadratic, dense$quadratic, tolerance = 1e-10)
  expect_equal(blockwise$logdet_v, dense$logdet_v, tolerance = 1e-10)
})

test_that("Satterthwaite derivatives are cached once per fit", {
  fit <- fit_orthodont_intercept()
  expect_length(fit$beta_vcov_derivatives, length(fit$eta))
  expect_true(all(vapply(
    fit$beta_vcov_derivatives,
    function(x) identical(dim(x), dim(vcov(fit, adjusted = FALSE))),
    logical(1L)
  )))

  cached <- fit
  cached$design <- NULL
  expect_equal(
    satterthwaite_df(cached, c(0, 1, 0)),
    satterthwaite_df(fit, c(0, 1, 0)),
    tolerance = 1e-12
  )
})

test_that("simulate.lmm is reproducible without changing the RNG state", {
  fit <- fit_orthodont_intercept(ddf = "residual")
  set.seed(2027)
  initial_seed <- .Random.seed
  first <- simulate(fit, nsim = 3, seed = 91)
  second <- simulate(fit, nsim = 3, seed = 91)

  expect_identical(.Random.seed, initial_seed)
  expect_equal(first, second)
  expect_s3_class(first, "data.frame")
  expect_equal(dim(first), c(nobs(fit), 3L))
  expect_named(first, paste0("sim_", 1:3))
  expect_error(simulate(fit, nsim = 0), "positive integer")
})

test_that("update.lmm changes model components while retaining the data", {
  fit <- fit_orthodont_intercept(ddf = "residual")
  updated <- update(
    fit,
    . ~ . - Sex,
    method = "ML",
    ddf = "residual"
  )

  expect_identical(updated$method, "ML")
  expect_equal(nobs(updated), nobs(fit))
  expect_named(fixef(updated), c("(Intercept)", "age"))
  expect_true(is.language(update(fit, evaluate = FALSE)))
})

test_that("predict.lmm supplies standard errors and intervals", {
  fit <- fit_orthodont_intercept()
  standard_errors <- predict(fit, re.form = NA, se.fit = TRUE)
  confidence <- predict(fit, re.form = NA, interval = "confidence")
  prediction <- predict(fit, re.form = NA, interval = "prediction")

  expect_named(standard_errors, c("fit", "se.fit", "df"))
  expect_true(all(is.finite(standard_errors$se.fit)))
  expect_identical(colnames(confidence), c("fit", "lwr", "upr"))
  expect_true(all(
    prediction[, "upr"] - prediction[, "lwr"] >
      confidence[, "upr"] - confidence[, "lwr"]
  ))
  expect_equal(predict(fit, re.form = ~0), fitted(fit, type = "marginal"))
  expect_equal(predict(fit, re.form = ~1 | Subject), fitted(fit))
  expect_error(
    predict(fit, re.form = ~age | Subject),
    "Unknown random-effect"
  )
})

test_that("known precision weights agree with weighted least squares", {
  data <- data.frame(
    x = 0:9,
    weight = seq(0.5, 2, length.out = 10),
    y = c(1.1, 1.8, 2.7, 3.0, 3.9, 4.8, 5.1, 6.2, 6.8, 7.7)
  )
  fit <- lmm(
    data,
    y ~ x,
    weights = weight,
    method = "ML",
    ddf = "residual"
  )
  reference <- stats::lm(y ~ x, data = data, weights = weight)

  expect_equal(fixef(fit), stats::coef(reference), tolerance = 1e-8)
  expect_equal(
    diag(fit$covariance$r) * data$weight,
    rep(VarCorr(fit)$estimate[[1L]], nrow(data)),
    tolerance = 1e-10
  )
  expect_error(
    lmm(data, y ~ x, weights = c(rep(1, 9), 0)),
    "positive finite"
  )
})

test_that("explicit offsets and contrasts are honored", {
  data <- data.frame(
    x = 0:11,
    group = factor(rep(c("a", "b", "c"), each = 4)),
    exposure = seq(0.2, 1.3, length.out = 12)
  )
  data$y <- 2 + 0.6 * data$x + data$exposure +
    rep(c(-0.2, 0.1, 0.3, -0.1), 3)
  fit <- lmm(
    data,
    y ~ x,
    offset = exposure,
    method = "ML",
    ddf = "residual"
  )
  reference <- stats::lm(y ~ x + offset(exposure), data = data)
  newdata <- data.frame(x = c(1, 2), exposure = c(0.5, 1.5))

  expect_equal(fixef(fit), stats::coef(reference), tolerance = 1e-8)
  expect_equal(
    predict(fit, newdata = newdata, re.form = NA),
    stats::predict(reference, newdata = newdata),
    tolerance = 1e-8
  )

  factor_fit <- lmm(
    data,
    y ~ group,
    contrasts = list(group = "contr.sum"),
    method = "ML",
    ddf = "residual"
  )
  expect_identical(
    colnames(model.matrix(factor_fit)),
    colnames(stats::model.matrix(
      y ~ group,
      data,
      contrasts.arg = list(group = "contr.sum")
    ))
  )
})

test_that("the multi-df fallback is explicit", {
  expect_warning(
    expect_identical(f_statistic_df(c(1.5, 5)), 2),
    "conservative fallback"
  )
})
