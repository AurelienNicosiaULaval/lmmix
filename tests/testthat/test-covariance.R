test_that("every residual covariance structure is positive definite", {
  for (structure in c("id", "cs", "ar1", "toep", "un")) {
    repeated <- if (structure == "id") NULL else ~ Time | Center:Subject
    fit <- lmm(
      multicentre,
      Y ~ Drug * Time,
      repeated = repeated,
      structure = structure
    )

    eigenvalues <- eigen(
      fit$covariance$v,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    expect_true(all(eigenvalues > 0), info = structure)
    expect_identical(fit$convergence$code, 0L, info = structure)
  }
})

test_that("AR1 accepts negative correlations", {
  data <- data.frame(
    subject = factor(rep(seq_len(20), each = 4)),
    time = rep(seq_len(4), 20)
  )
  set.seed(42)
  innovation <- stats::rnorm(nrow(data))
  response <- numeric(nrow(data))
  for (subject in levels(data$subject)) {
    rows <- which(data$subject == subject)
    response[rows[[1L]]] <- innovation[rows[[1L]]] / sqrt(1 - 0.5^2)
    for (index in 2:length(rows)) {
      response[rows[[index]]] <-
        -0.5 * response[rows[[index - 1L]]] + innovation[rows[[index]]]
    }
  }
  data$response <- response
  fit <- lmm(
    data,
    response ~ 1,
    repeated = ~ time | subject,
    structure = "ar1"
  )

  correlation <- VarCorr(fit)$estimate[VarCorr(fit)$component == "cor"]
  expect_lt(correlation, 0)
})

test_that("unstructured covariance has occasion-specific variances", {
  fit <- lmm(
    multicentre,
    Y ~ Drug * Time,
    repeated = ~ Time | Center:Subject,
    structure = "un"
  )
  covariance <- VarCorr(fit)

  expect_equal(sum(covariance$component == "var"), 3L)
  expect_equal(sum(covariance$component == "cor"), 3L)
})
