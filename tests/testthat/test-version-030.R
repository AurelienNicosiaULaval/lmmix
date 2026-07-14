test_that("fixed-band Toeplitz reproduces SAS Example 79.6", {
  fit <- lmm(
    sas_line_source,
    Y ~ (Cult + Dir + Irrig)^2,
    random = list(
      Block = ~ 1 | Block,
      Block.Dir = ~ 1 | Block:Dir,
      Block.Irrig = ~ 1 | Block:Irrig
    ),
    repeated = ~ Sbplt | Block:Cult,
    structure = "toep(4)",
    ddf = "residual"
  )

  expected <- c(
    0.2194,
    0.01768,
    0.03539,
    0.2850,
    0.007986 / 0.2850,
    0.001452 / 0.2850,
    -0.09253 / 0.2850
  )
  expect_identical(fit$structure, "toep")
  expect_identical(fit$structure_label, "toep(4)")
  expect_identical(fit$covariance_order, 4L)
  expect_equal(VarCorr(fit)$estimate, expected, tolerance = 5e-4)
  expect_equal(deviance(fit), 183.23797748, tolerance = 1e-5)
  expect_equal(fit$covariance$residual_base[1, 5], 0, tolerance = 1e-12)
  expect_message(print(fit), "TOEP\\(4\\)")
})

test_that("fixed-band Toeplitz validates its order", {
  expect_error(
    lmm(
      sas_growth,
      y ~ Gender * Age,
      repeated = ~ Age | Person,
      structure = "toep(5)"
    ),
    "cannot exceed"
  )
  expect_error(
    lmm(
      sas_growth,
      y ~ Gender * Age,
      repeated = ~ Age | Person,
      structure = "toep(0)"
    ),
    "must be one of"
  )
})

test_that("covariance confidence intervals respect parameter bounds", {
  fit <- fit_orthodont_intercept()
  covariance_intervals <- confint(fit, parm = "theta_")
  all_intervals <- confint(fit, parm = c("beta_", "theta_"))

  expect_equal(nrow(covariance_intervals), nrow(VarCorr(fit)))
  expect_true(all(covariance_intervals[, "Lower"] > 0))
  expect_true(all(covariance_intervals[, "Upper"] > covariance_intervals[, 1]))
  expect_equal(
    nrow(all_intervals),
    length(fixef(fit)) + nrow(VarCorr(fit))
  )

  slope_fit <- lmm(
    orthodont_data(),
    distance ~ age + Sex,
    random = ~ age | Subject
  )
  correlation_name <- grep(
    "random.cor",
    rownames(confint(slope_fit, parm = "theta_")),
    value = TRUE
  )
  correlation_interval <- confint(slope_fit, parm = correlation_name)
  expect_true(all(correlation_interval > -1 & correlation_interval < 1))
  expect_error(confint(fit, parm = "unknown"), "Unknown parameter")
  expect_error(confint(fit, level = 1), "strictly between")
})

test_that("parametric-bootstrap likelihood-ratio tests are reproducible", {
  data <- orthodont_data()
  marginal <- lmm(data, distance ~ age + Sex, method = "ML")
  random <- lmm(
    data,
    distance ~ age + Sex,
    random = ~ 1 | Subject,
    method = "ML"
  )

  set.seed(2026)
  initial_seed <- .Random.seed
  first <- anova(
    marginal,
    random,
    test = "parametric.bootstrap",
    nsim = 5,
    seed = 123
  )
  second <- anova(
    marginal,
    random,
    test = "parametric.bootstrap",
    nsim = 5,
    seed = 123
  )

  expect_identical(.Random.seed, initial_seed)
  expect_equal(first$p.value, second$p.value)
  expect_identical(attr(first, "test"), "parametric.bootstrap")
  expect_identical(attr(first, "nsim"), 5L)
  expect_identical(attr(first, "bootstrap")[[2L]]$successful, 5L)
  expect_true(first$p.value[[2L]] %in% (seq_len(6L) / 6))
  expect_message(print(first), "bootstrap simulations")
  expect_error(
    anova(marginal, random, test = "parametric.bootstrap", nsim = 0),
    "positive integer"
  )
  expect_error(
    anova(
      marginal,
      random,
      test = "parametric.bootstrap",
      nsim = 2,
      seed = -1
    ),
    "non-negative integer"
  )
  expect_error(
    anova(marginal, test = "parametric.bootstrap"),
    "two or more models"
  )
})
