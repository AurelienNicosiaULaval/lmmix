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

  expect_absolute_error_below(
    unname(result$estimate),
    unname(reference_table[, "Estimate"]),
    tolerance = 1e-5
  )
  expect_absolute_error_below(
    unname(result$std.error),
    unname(reference_table[, "Std. Error"]),
    tolerance = 1e-5
  )
  expect_absolute_error_below(
    unname(result$df),
    unname(reference_table[, "df"]),
    tolerance = 5e-3
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

  expect_absolute_error_below(
    unname(result$num.df),
    unname(reference_table[, "NumDF"]),
    tolerance = 1e-6
  )
  expect_absolute_error_below(
    unname(result$den.df),
    unname(reference_table[, "DenDF"]),
    tolerance = 5e-3
  )
  expect_absolute_error_below(
    unname(result$statistic),
    unname(reference_table[, "F value"]),
    tolerance = 5e-3
  )
})

test_that("residual degrees of freedom are available", {
  fit <- fit_orthodont_intercept(ddf = "residual")
  result <- generics::tidy(fit)

  expect_equal(unique(result$df), nobs(fit) - length(fixef(fit)))
})

test_that("the ANOVA print method reports its denominator df method", {
  fit <- fit_orthodont_intercept()
  expect_output(print(anova(fit)), "satterthwaite")
})

test_that("Kenward-Roger agrees with lmerTest degrees of freedom", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("pbkrtest")
  fit <- fit_orthodont_intercept(ddf = "kenward-roger")
  reference <- lmerTest::lmer(
    distance ~ age + Sex + (1 | Subject),
    data = orthodont_data(),
    REML = TRUE
  )
  reference_df <- vapply(seq_along(fixef(fit)), function(index) {
    contrast <- numeric(length(fixef(fit)))
    contrast[[index]] <- 1
    lmerTest::contest1D(
      reference,
      contrast,
      ddf = "Kenward-Roger"
    )$df
  }, numeric(1L))

  expect_equal(summary(fit)$fixed$df, reference_df, tolerance = 5e-3)
  expect_true(all(eigen(vcov(fit), symmetric = TRUE)$values > 0))
  expect_false(isTRUE(all.equal(vcov(fit), vcov(fit, adjusted = FALSE))))
})

test_that("Kenward-Roger agrees with mmrm for marginal AR1 models", {
  skip_if_not_installed("mmrm")
  data <- orthodont_data()
  fit <- lmm(
    data,
    distance ~ age + Sex,
    repeated = ~ Occasion | Subject,
    structure = "ar1",
    ddf = "kenward-roger"
  )
  reference_formula <- stats::as.formula(
    "distance ~ age + Sex + ar1(Occasion | Subject)",
    env = asNamespace("mmrm")
  )
  reference <- mmrm::mmrm(
    reference_formula,
    data = data,
    method = "Kenward-Roger"
  )
  reference_table <- coef(summary(reference))

  expect_equal(fixef(fit), coef(reference), tolerance = 1e-5)
  expect_equal(
    summary(fit)$fixed$std.error,
    unname(reference_table[, "Std. Error"]),
    tolerance = 3e-3
  )
  expect_equal(
    summary(fit)$fixed$df,
    unname(reference_table[, "df"]),
    tolerance = 5e-3
  )
})

test_that("Kenward-Roger requires REML", {
  expect_error(
    lmm(
      orthodont_data(),
      distance ~ age + Sex,
      random = ~ 1 | Subject,
      method = "ML",
      ddf = "kenward-roger"
    ),
    "requires.*REML"
  )
})
