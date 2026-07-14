crossed_random_data <- function() {
  set.seed(42)
  data <- expand.grid(
    site = factor(seq_len(8)),
    observer = factor(seq_len(6)),
    replicate = seq_len(3)
  )
  site_effect <- stats::rnorm(8, sd = 0.8)
  observer_effect <- stats::rnorm(6, sd = 0.5)
  data$response <- 2 +
    site_effect[data$site] +
    observer_effect[data$observer] +
    stats::rnorm(nrow(data), sd = 0.4)
  data
}

test_that("multiple crossed random terms agree with lme4", {
  skip_if_not_installed("lme4")
  data <- crossed_random_data()
  fit <- lmm(
    data,
    response ~ 1,
    random = list(site = ~ 1 | site, observer = ~ 1 | observer)
  )
  reference <- lme4::lmer(
    response ~ 1 + (1 | site) + (1 | observer),
    data = data,
    REML = TRUE
  )
  reference_variance <- as.data.frame(lme4::VarCorr(reference))$vcov

  expect_equal(
    as.numeric(logLik(fit)),
    as.numeric(logLik(reference)),
    tolerance = 1e-5
  )
  expect_equal(VarCorr(fit)$estimate, reference_variance, tolerance = 1e-4)
  expect_named(ranef(fit), c("site", "observer"))
  expect_true(all(vapply(ranef(fit), inherits, logical(1L), what = "tbl_df")))
})

test_that("Kenward-Roger supports multiple random terms", {
  data <- crossed_random_data()
  fit <- lmm(
    data,
    response ~ 1,
    random = list(site = ~ 1 | site, observer = ~ 1 | observer),
    ddf = "kenward-roger"
  )
  fixed <- generics::tidy(fit)

  expect_true(all(is.finite(fixed$std.error)))
  expect_true(all(is.finite(fixed$df)))
  expect_true(all(eigen(vcov(fit), symmetric = TRUE)$values > 0))
})

test_that("random formula lists are validated", {
  data <- crossed_random_data()
  expect_error(lmm(data, response ~ 1, random = list(1)), "Every element")
  expect_error(
    lmm(data, response ~ 1, random = list(~ 1 | site, ~ 1 | site)),
    "Duplicated"
  )
})

test_that("fixed-effect model comparisons refit REML models with ML", {
  skip_if_not_installed("lme4")
  data <- orthodont_data()
  reduced <- lmm(data, distance ~ age, random = ~ 1 | Subject)
  full <- lmm(data, distance ~ age + Sex, random = ~ 1 | Subject)
  expect_message(result <- anova(reduced, full), "Refitting")

  reference_reduced <- lme4::lmer(
    distance ~ age + (1 | Subject),
    data = data,
    REML = FALSE
  )
  reference_full <- lme4::lmer(
    distance ~ age + Sex + (1 | Subject),
    data = data,
    REML = FALSE
  )
  reference <- anova(reference_reduced, reference_full, refit = FALSE)

  expect_s3_class(result, "anova.lmm_list")
  expect_true(isTRUE(attr(result, "refitted")))
  expect_equal(result$Chisq[[2L]], reference$Chisq[[2L]], tolerance = 1e-5)
  expect_equal(
    result$p.value[[2L]],
    reference$`Pr(>Chisq)`[[2L]],
    tolerance = 1e-6
  )
  expect_error(anova(reduced, full, refit = FALSE), "must be compared with ML")
})

test_that("covariance comparisons require ML and warn about boundaries", {
  data <- orthodont_data()
  marginal <- lmm(data, distance ~ age + Sex, method = "ML")
  random <- lmm(
    data,
    distance ~ age + Sex,
    random = ~ 1 | Subject,
    method = "ML"
  )
  expect_warning(result <- anova(marginal, random), "boundary")
  expect_true(is.finite(result$Chisq[[2L]]))

  marginal_reml <- lmm(data, distance ~ age + Sex)
  random_reml <- lmm(data, distance ~ age + Sex, random = ~ 1 | Subject)
  expect_error(anova(marginal_reml, random_reml), "fitted with ML")
})

test_that("na.exclude restores omitted rows", {
  data <- orthodont_data()
  data$distance[c(2L, 5L)] <- NA_real_
  fit <- lmm(
    data,
    distance ~ age + Sex,
    random = ~ 1 | Subject,
    na.action = stats::na.exclude
  )

  expect_length(fitted(fit), nrow(data))
  expect_length(residuals(fit), nrow(data))
  expect_true(all(is.na(fitted(fit)[c(2L, 5L)])))
  expect_equal(nrow(generics::augment(fit)), nrow(data))
  expect_error(
    lmm(
      data,
      distance ~ age + Sex,
      random = ~ 1 | Subject,
      na.action = stats::na.fail
    ),
    "Missing values"
  )
})

test_that("adjusted pairwise results include simultaneous intervals", {
  skip_if_not_installed("emmeans")
  fit <- fit_multicentre_ar1()
  pointwise <- lsmeans(fit, pairwise ~ Drug)$contrasts
  adjusted <- lsmeans(
    fit,
    pairwise ~ Drug,
    adjust = "holm",
    conf_adjust = "auto"
  )$contrasts

  pointwise_width <- pointwise$conf.high - pointwise$conf.low
  adjusted_width <- adjusted$conf.high - adjusted$conf.low
  expect_true(all(adjusted_width > pointwise_width))
  expect_identical(attr(adjusted, "p.adjust"), "holm")
  expect_identical(attr(adjusted, "conf.adjust"), "bonferroni")

  bonferroni <- lsmeans(
    fit,
    pairwise ~ Drug,
    adjust = "bonferroni"
  )$contrasts
  reference <- suppressMessages(emmeans::emmeans(fit, ~Drug)) |>
    emmeans::contrast(method = "pairwise", adjust = "bonferroni") |>
    confint() |>
    as.data.frame()
  expect_equal(bonferroni$conf.low, reference$lower.CL, tolerance = 1e-6)
  expect_equal(bonferroni$conf.high, reference$upper.CL, tolerance = 1e-6)
})

test_that("automatic optimization records attempts", {
  fit <- fit_orthodont_intercept()
  expect_s3_class(fit$convergence$attempts, "tbl_df")
  expect_true(any(fit$convergence$attempts$selected))
  expect_error(lmm_control(max_restarts = -1), "non-negative")
  expect_error(lmm_control(restart_scale = -1), "non-negative")

  restarted <- suppressWarnings(lmm(
    orthodont_data(),
    distance ~ age + Sex,
    random = ~ 1 | Subject,
    control = lmm_control(max_iter = 1, max_restarts = 2)
  ))
  expect_equal(nrow(restarted$convergence$attempts), 3L)
  expect_identical(
    restarted$convergence$attempts$optimizer,
    c("nlminb", "optim", "optim")
  )
})
