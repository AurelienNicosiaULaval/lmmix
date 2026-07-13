test_that("core S3 methods return documented shapes", {
  fit <- fit_orthodont_intercept()

  expect_message(print(fit), "Linear mixed model")
  model_summary <- summary(fit)
  expect_s3_class(model_summary, "summary.lmm")
  expect_message(print(model_summary), "Fixed effects")
  expect_equal(coef(fit), fixef(fit))
  expect_named(fixef(fit), c("(Intercept)", "age", "SexFemale"))
  expect_s3_class(VarCorr(fit), "tbl_df")
  expect_s3_class(ranef(fit), "tbl_df")
  expect_equal(dim(vcov(fit)), c(3L, 3L))
  expect_s3_class(logLik(fit), "logLik")
  expect_length(fitted(fit), nobs(fit))
  expect_length(fitted(fit, type = "marginal"), nobs(fit))
  expect_length(residuals(fit), nobs(fit))
  expect_length(residuals(fit, type = "pearson"), nobs(fit))
  expect_length(residuals(fit, type = "marginal"), nobs(fit))
  expect_length(predict(fit), nobs(fit))
  expect_equal(predict(fit, re.form = NA), fitted(fit, type = "marginal"))
  expect_equal(model.matrix(fit), fit$design$x)
  expect_equal(model.frame(fit), fit$model_frame)
  expect_identical(formula(fit), fit$formula)
  expect_named(
    formula(fit, fixed.only = FALSE),
    c("fixed", "random", "repeated")
  )
  expect_identical(terms(fit), fit$terms)
  type3 <- anova(fit)
  expect_s3_class(type3, "anova.lmm")
  expect_message(print(type3), "Type III")
  expect_equal(dim(confint(fit)), c(3L, 2L))
  expect_equal(dim(confint(fit, parm = 2)), c(1L, 2L))
  expect_error(confint(fit, parm = "unknown"), "Unknown fixed-effect")
  expect_error(anova(fit, type = 2), "Only type III")
  expect_error(anova(fit, fit), "Model-comparison")
})

test_that("broom methods return tibbles", {
  fit <- fit_orthodont_intercept()
  tidy_result <- generics::tidy(fit)
  glance_result <- generics::glance(fit)
  augment_result <- generics::augment(fit)

  expect_s3_class(tidy_result, "tbl_df")
  expect_named(
    tidy_result,
    c("effect", "term", "estimate", "std.error", "statistic", "df", "p.value")
  )
  expect_s3_class(glance_result, "tbl_df")
  expect_equal(nrow(glance_result), 1L)
  expect_named(
    glance_result,
    c("logLik", "AIC", "BIC", "deviance", "df", "nobs", "convergence", "method")
  )
  expect_s3_class(augment_result, "tbl_df")
  expected_columns <- c(".fitted", ".resid", ".std.resid")
  expect_true(all(expected_columns %in% names(augment_result)))

  covariance_tidy <- generics::tidy(fit, effects = "ran_pars")
  all_tidy <- generics::tidy(fit, effects = "all", conf.int = TRUE)
  expect_true(all(covariance_tidy$effect == "ran_pars"))
  expect_true(all(c("conf.low", "conf.high") %in% names(all_tidy)))
})

test_that("prediction handles known and new random-effect groups", {
  fit <- fit_orthodont_intercept()
  newdata <- data.frame(
    age = c(10, 10),
    Sex = factor(c("Male", "Male"), levels = c("Male", "Female")),
    Subject = factor(c("M01", "NEW"))
  )
  conditional <- predict(fit, newdata = newdata)
  marginal <- predict(fit, newdata = newdata, re.form = NA)

  expect_false(isTRUE(all.equal(conditional[[1L]], marginal[[1L]])))
  expect_equal(conditional[[2L]], marginal[[2L]])
  expect_equal(nrow(model.matrix(fit, data = newdata)), 2L)

  augmented <- generics::augment(
    fit,
    newdata = newdata[c("age", "Sex", "Subject")]
  )
  expect_true(all(is.na(augmented$.resid)))
})

test_that("missing responses are removed reproducibly", {
  fit <- fit_multicentre_ar1()
  augmented <- generics::augment(fit)

  expect_equal(nrow(augmented), sum(!is.na(multicentre$Y)))
  expect_false(anyNA(augmented$Y))
})

test_that("marginal models have no random-effect table", {
  fit <- lmm(
    orthodont_data(),
    distance ~ age + Sex,
    repeated = ~ Occasion | Subject,
    structure = "ar1"
  )

  expect_equal(nrow(ranef(fit)), 0L)
})
