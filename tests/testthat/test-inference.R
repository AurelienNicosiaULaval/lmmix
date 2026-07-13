test_that("one-dimensional Satterthwaite inference agrees with lmerTest", {
  skip_if_not_installed("lmerTest")
  data <- orthodont_data()
  fit <- fit_orthodont_intercept()
  reference <- lmerTest::lmer(
    distance ~ age + Sex + (1 | Subject),
    data = data,
    REML = TRUE
  )
  reference_table <- summary(reference)$coefficients
  result <- generics::tidy(fit)

  expect_equal(
    unname(result$estimate),
    unname(reference_table[, "Estimate"]),
    tolerance = 1e-5
  )
  expect_equal(
    unname(result$std.error),
    unname(reference_table[, "Std. Error"]),
    tolerance = 1e-5
  )
  expect_equal(
    unname(result$df),
    unname(reference_table[, "df"]),
    tolerance = 1e-3
  )
})

test_that("type III tests agree with lmerTest", {
  skip_if_not_installed("lmerTest")
  data <- orthodont_data()
  fit <- fit_orthodont_intercept()
  reference <- lmerTest::lmer(
    distance ~ age + Sex + (1 | Subject),
    data = data,
    REML = TRUE
  )
  reference_table <- anova(reference, type = 3)
  result <- anova(fit)

  expect_equal(
    unname(result$num.df),
    unname(reference_table[, "NumDF"]),
    tolerance = 1e-6
  )
  expect_equal(
    unname(result$den.df),
    unname(reference_table[, "DenDF"]),
    tolerance = 1e-3
  )
  expect_equal(
    unname(result$statistic),
    unname(reference_table[, "F value"]),
    tolerance = 1e-3
  )
})

test_that("residual degrees of freedom are available", {
  fit <- fit_orthodont_intercept(ddf = "residual")
  result <- generics::tidy(fit)

  expect_equal(unique(result$df), nobs(fit) - length(fixef(fit)))
})

test_that("Kenward-Roger is rejected rather than silently approximated", {
  expect_error(
    fit_orthodont_intercept(ddf = "kenward-roger"),
    "not implemented"
  )
})
