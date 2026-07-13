test_that("LS-means average nuisance-factor levels equally", {
  fit <- fit_multicentre_ar1()
  result <- lsmeans(fit, ~Drug)

  expect_s3_class(result, "lmm_lsmeans")
  expect_equal(nrow(result), 3L)
  expect_equal(
    result$estimate,
    c(13.973956, 19.710671, 16.900468),
    tolerance = 1e-3
  )
})

test_that("pairwise LS-mean contrasts are automatic", {
  fit <- fit_multicentre_ar1()
  result <- lsmeans(fit, pairwise ~ Drug)

  expect_s3_class(result, "lmm_lsmeans_list")
  expect_equal(nrow(result$contrasts), 3L)
  expect_equal(result$contrasts$contrast, c("1 - 2", "1 - 3", "2 - 3"))
  expect_equal(
    result$contrasts$estimate,
    c(-5.736715, -2.926512, 2.810203),
    tolerance = 1e-3
  )
  expect_message(print(result), "Pairwise contrasts")
})

test_that("emmeans uses the lmm basis and Satterthwaite df", {
  skip_if_not_installed("emmeans")
  fit <- fit_orthodont_intercept()
  internal <- lsmeans(fit, ~Sex)
  external <- as.data.frame(emmeans::emmeans(fit, ~Sex))

  expect_equal(internal$estimate, external$emmean, tolerance = 1e-6)
  expect_equal(internal$std.error, external$SE, tolerance = 1e-6)
  expect_equal(internal$df, external$df, tolerance = 1e-3)
})

test_that("character specs and reference values are supported", {
  fit <- fit_orthodont_intercept()
  result <- lsmeans(fit, "Sex", at = list(age = 12))

  expect_equal(nrow(result), 2L)
  expect_equal(result$estimate, c(25.62894, 23.30792), tolerance = 1e-3)
})

test_that("LS-mean inputs are validated", {
  fit <- fit_orthodont_intercept()

  expect_error(lsmeans(fit, ~unknown), "Unknown reference-grid")
  expect_error(lsmeans(fit, ~Sex, at = 1), "named list")
  expect_error(lsmeans(fit, ~Sex, level = 1), "strictly between")
  expect_error(lsmeans(fit, ~Sex, adjust = "invalid"), "Unknown")
  expect_error(lsmeans(fit, response ~ Sex), "left side")
  expect_error(lsmeans(fit, character()), "at least one")
})
