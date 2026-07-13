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

  expect_absolute_error_below(
    fixef(fit),
    nlme::fixef(reference),
    tolerance = 1e-6
  )
  expect_absolute_error_below(
    as.numeric(logLik(fit)),
    as.numeric(logLik(reference)),
    tolerance = 1e-6
  )
  reference_covariance <- as.matrix(nlme::getVarCov(
    reference,
    type = "random.effects"
  ))
  expected_covariance <- c(
    reference_covariance[1L, 1L],
    reference$sigma^2
  )
  expect_absolute_error_below(
    VarCorr(fit)$estimate,
    expected_covariance,
    tolerance = 1e-4
  )
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

  expect_absolute_error_below(
    fixef(fit),
    nlme::fixef(reference),
    tolerance = 1e-5
  )
  expect_absolute_error_below(
    as.numeric(logLik(fit)),
    as.numeric(logLik(reference)),
    tolerance = 1e-5
  )
  reference_covariance <- as.matrix(nlme::getVarCov(
    reference,
    type = "random.effects"
  ))
  expected_covariance <- c(
    reference_covariance[1L, 1L],
    reference_covariance[2L, 2L],
    reference_covariance[1L, 2L] / sqrt(
      reference_covariance[1L, 1L] * reference_covariance[2L, 2L]
    ),
    reference$sigma^2
  )
  expect_absolute_error_below(
    VarCorr(fit)$estimate,
    expected_covariance,
    tolerance = 1e-4
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

  expect_absolute_error_below(
    fixef(fit),
    nlme::fixef(reference),
    tolerance = 1e-6
  )
  expect_absolute_error_below(
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

    expect_absolute_error_below(
      fixef(fit),
      stats::coef(reference),
      tolerance = 1e-5
    )
    expect_absolute_error_below(
      as.numeric(logLik(fit)),
      as.numeric(logLik(reference)),
      tolerance = 1e-5
    )
    expected_covariance <- c(
      residual.var = reference$sigma^2,
      residual.cor = as.numeric(
        stats::coef(reference$modelStruct$corStruct, unconstrained = FALSE)
      )
    )
    expect_absolute_error_below(
      stats::setNames(VarCorr(fit)$estimate, c("residual.var", "residual.cor")),
      expected_covariance,
      tolerance = 1e-5
    )
  }
})

test_that("marginal AR1 fit agrees with mmrm", {
  skip_if_not_installed("mmrm")
  data <- orthodont_data()
  fit <- lmm(
    data,
    distance ~ age + Sex,
    repeated = ~ Occasion | Subject,
    structure = "ar1"
  )
  reference <- mmrm::mmrm(
    distance ~ age + Sex,
    data = data,
    covariance = mmrm::cov_struct(
      type = "ar1",
      visits = "Occasion",
      subject = "Subject"
    ),
    reml = TRUE
  )

  expect_absolute_error_below(
    fixef(fit),
    stats::coef(reference),
    tolerance = 1e-5
  )
  expect_absolute_error_below(
    as.numeric(logLik(fit)),
    as.numeric(logLik(reference)),
    tolerance = 1e-5
  )
  expect_absolute_error_below(
    unname(fit$covariance$residual_base),
    unname(mmrm::component(reference, "varcor")),
    tolerance = 1e-4
  )
})

test_that("the combined random-effect and AR1 model converges", {
  fit <- fit_multicentre_ar1()

  expect_s3_class(fit, "lmm")
  expect_identical(fit$convergence$code, 0L)
  expect_true(fit$convergence$hessian_positive_definite)
  expect_true(all(eigen(fit$covariance$v, symmetric = TRUE)$values > 0))
  expect_equal(nobs(fit), sum(!is.na(multicentre$Y)))
})

test_that("every residual structure fits in a combined model", {
  for (structure in c("id", "cs", "ar1", "toep", "un")) {
    expect_no_warning(
      fit <- lmm(
        multicentre,
        Y ~ Drug * Time,
        random = ~ 1 | Center,
        repeated = ~ Time | Center:Drug:Subject,
        structure = structure
      )
    )

    expect_identical(fit$convergence$code, 0L)
    expect_true(fit$convergence$hessian_positive_definite)
    expect_gt(
      min(eigen(
        fit$covariance$v,
        symmetric = TRUE,
        only.values = TRUE
      )$values),
      0
    )
  }
})
