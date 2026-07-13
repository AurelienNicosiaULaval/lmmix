contrast_variance <- function(object, contrast, eta = object$eta) {
  beta_vcov <- beta_vcov_at_eta(eta, object$design)
  drop(contrast %*% beta_vcov %*% contrast)
}

satterthwaite_df <- function(object, contrast) {
  contrast <- as.numeric(contrast)
  variance <- contrast_variance(object, contrast)
  invalid_variance <- !is.finite(variance) || variance <= 0
  if (invalid_variance || any(!is.finite(object$eta_vcov))) {
    return(NA_real_)
  }

  gradient <- numDeriv::grad(
    func = function(eta) contrast_variance(object, contrast, eta),
    x = object$eta,
    method = object$control$deriv_method
  )
  variance_of_variance <- drop(
    crossprod(gradient, object$eta_vcov %*% gradient)
  )
  if (!is.finite(variance_of_variance) || variance_of_variance <= 0) {
    return(Inf)
  }

  2 * variance^2 / variance_of_variance
}

contrast_df <- function(object, contrast) {
  if (object$ddf == "residual") {
    return(nrow(object$design$x) - ncol(object$design$x))
  }
  satterthwaite_df(object, contrast)
}

contrast_statistics <- function(object, contrast, level = 0.95) {
  contrast <- as.numeric(contrast)
  estimate <- drop(contrast %*% object$coefficients)
  variance <- drop(contrast %*% object$beta_vcov %*% contrast)
  standard_error <- sqrt(max(variance, 0))
  df <- contrast_df(object, contrast)
  statistic <- estimate / standard_error
  p_value <- 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE)
  critical <- stats::qt((1 + level) / 2, df = df)

  c(
    estimate = estimate,
    std.error = standard_error,
    statistic = statistic,
    df = df,
    p.value = p_value,
    conf.low = estimate - critical * standard_error,
    conf.high = estimate + critical * standard_error
  )
}

fixed_effects_table <- function(object, level = 0.95) {
  p <- length(object$coefficients)
  statistics <- lapply(seq_len(p), function(index) {
    contrast <- numeric(p)
    contrast[[index]] <- 1
    contrast_statistics(object, contrast, level = level)
  })
  statistics <- do.call(rbind, statistics)

  tibble::tibble(
    term = names(object$coefficients),
    estimate = unname(statistics[, "estimate"]),
    std.error = unname(statistics[, "std.error"]),
    statistic = unname(statistics[, "statistic"]),
    df = unname(statistics[, "df"]),
    p.value = unname(statistics[, "p.value"]),
    conf.low = unname(statistics[, "conf.low"]),
    conf.high = unname(statistics[, "conf.high"])
  )
}

sum_contrast_matrix <- function(object) {
  factor_columns <- vapply(object$model_frame, is.factor, logical(1L))
  factors <- object$model_frame[factor_columns]
  contrasts_sum <- lapply(factors, function(x) {
    if (nlevels(x) > 1L) stats::contr.sum(nlevels(x)) else NULL
  })
  contrasts_sum <- contrasts_sum[!vapply(contrasts_sum, is.null, logical(1L))]

  x_sum <- stats::model.matrix(
    object$terms,
    data = object$model_frame,
    contrasts.arg = contrasts_sum
  )
  transformation <- qr.solve(object$design$x, x_sum, tol = 1e-10)
  inverse_transformation <- solve(transformation)

  list(
    matrix = x_sum,
    inverse_transformation = inverse_transformation,
    assignment = attr(x_sum, "assign")
  )
}

type3_contrasts <- function(object) {
  sum_matrix <- sum_contrast_matrix(object)
  term_labels <- attr(object$terms, "term.labels")

  stats::setNames(
    lapply(seq_along(term_labels), function(index) {
      rows <- which(sum_matrix$assignment == index)
      sum_matrix$inverse_transformation[rows, , drop = FALSE]
    }),
    term_labels
  )
}

f_statistic_df <- function(component_df, tolerance = 1e-8) {
  if (length(component_df) == 1L) {
    return(component_df)
  }
  if (all(abs(diff(component_df)) < tolerance)) {
    return(mean(component_df))
  }
  if (any(component_df <= 2)) {
    return(2)
  }

  expectation <- sum(component_df / (component_df - 2))
  2 * expectation / (expectation - length(component_df))
}

multi_df_test <- function(object, contrast) {
  covariance <- symmetrize(contrast %*% object$beta_vcov %*% t(contrast))
  eig <- eigen(covariance, symmetric = TRUE)
  tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(eig$values)))
  keep <- eig$values > tolerance
  numerator_df <- sum(keep)
  if (numerator_df == 0L) {
    return(c(num.df = 0, den.df = NA, statistic = NA, p.value = NA))
  }

  vectors <- eig$vectors[, keep, drop = FALSE]
  values <- eig$values[keep]
  projected_contrasts <- crossprod(vectors, contrast)
  projected_estimates <- drop(projected_contrasts %*% object$coefficients)
  statistic <- sum(projected_estimates^2 / values) / numerator_df

  if (object$ddf == "residual") {
    denominator_df <- nrow(object$design$x) - ncol(object$design$x)
  } else {
    component_df <- vapply(
      seq_len(numerator_df),
      function(index) {
        satterthwaite_df(object, projected_contrasts[index, ])
      },
      numeric(1L)
    )
    denominator_df <- f_statistic_df(component_df)
  }
  p_value <- stats::pf(
    statistic,
    df1 = numerator_df,
    df2 = denominator_df,
    lower.tail = FALSE
  )

  c(
    num.df = numerator_df,
    den.df = denominator_df,
    statistic = statistic,
    p.value = p_value
  )
}

type3_table <- function(object) {
  contrasts <- type3_contrasts(object)
  tests <- lapply(contrasts, function(contrast) {
    multi_df_test(object, contrast)
  })
  tests <- do.call(rbind, tests)

  tibble::tibble(
    term = names(contrasts),
    num.df = unname(tests[, "num.df"]),
    den.df = unname(tests[, "den.df"]),
    statistic = unname(tests[, "statistic"]),
    p.value = unname(tests[, "p.value"])
  )
}
