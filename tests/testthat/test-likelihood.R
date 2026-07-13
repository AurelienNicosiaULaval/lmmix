test_that("random-intercept REML agrees with nlme", {
  skip_if_not_installed("nlme")
  data <- orthodont_data()
  fit <- fit_orthodont_intercept()
  reference <- nlme::lme(
    distance ~ age + Sex,
    random = ~ 1 | Subject,
    data = data,
    method = "REML"
  )

  expect_equal(fixef(fit), nlme::fixef(reference), tolerance = 1e-6)
  expect_equal(
    as.numeric(logLik(fit)),
    as.numeric(logLik(reference)),
    tolerance = 1e-6
  )
  expect_equal(VarCorr(fit)$estimate, c(3.266784, 2.049456), tolerance = 1e-5)
})

test_that("random-slope REML agrees with nlme", {
  skip_if_not_installed("nlme")
  data <- orthodont_data()
  fit <- lmm(
    data,
    distance ~ age + Sex,
    random = ~ 1 + age | Subject
  )
  reference <- nlme::lme(
    distance ~ age + Sex,
    random = ~ age | Subject,
    data = data,
    method = "REML"
  )

  expect_equal(fixef(fit), nlme::fixef(reference), tolerance = 1e-5)
  expect_equal(
    as.numeric(logLik(fit)),
    as.numeric(logLik(reference)),
    tolerance = 1e-5
  )
})

test_that("ML agrees with nlme", {
  skip_if_not_installed("nlme")
  data <- orthodont_data()
  fit <- lmm(
    data,
    distance ~ age + Sex,
    random = ~ 1 | Subject,
    method = "ML"
  )
  reference <- nlme::lme(
    distance ~ age + Sex,
    random = ~ 1 | Subject,
    data = data,
    method = "ML"
  )

  expect_equal(fixef(fit), nlme::fixef(reference), tolerance = 1e-6)
  expect_equal(
    as.numeric(logLik(fit)),
    as.numeric(logLik(reference)),
    tolerance = 1e-6
  )
})

test_that("marginal CS and AR1 fits agree with nlme", {
  skip_if_not_installed("nlme")
  data <- orthodont_data()

  for (structure in c("cs", "ar1")) {
    fit <- lmm(
      data,
      distance ~ age + Sex,
      repeated = ~ Occasion | Subject,
      structure = structure
    )
    correlation <- if (structure == "cs") {
      nlme::corCompSymm(form = ~ 1 | Subject)
    } else {
      nlme::corAR1(form = ~ 1 | Subject)
    }
    reference <- nlme::gls(
      distance ~ age + Sex,
      data = data,
      correlation = correlation,
      method = "REML"
    )

    expect_equal(fixef(fit), stats::coef(reference), tolerance = 1e-5)
    expect_equal(
      as.numeric(logLik(fit)),
      as.numeric(logLik(reference)),
      tolerance = 1e-5
    )
  }
})

test_that("the combined random-effect and AR1 model converges", {
  fit <- fit_multicentre_ar1()

  expect_s3_class(fit, "lmm")
  expect_identical(fit$convergence$code, 0L)
  expect_true(fit$convergence$hessian_positive_definite)
  expect_true(all(eigen(fit$covariance$v, symmetric = TRUE)$values > 0))
  expect_equal(nobs(fit), sum(!is.na(multicentre$Y)))
})
