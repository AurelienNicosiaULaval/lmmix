covariance_derivatives <- function(eta, design) {
  n_parameters <- length(eta)
  base <- build_covariance(eta, design)$v
  step <- .Machine$double.eps^(1 / 4) * (abs(eta) + 1)
  first <- vector("list", n_parameters)
  second <- matrix(
    vector("list", n_parameters * n_parameters),
    nrow = n_parameters,
    ncol = n_parameters
  )
  plus <- vector("list", n_parameters)
  minus <- vector("list", n_parameters)

  for (index in seq_len(n_parameters)) {
    eta_plus <- eta
    eta_minus <- eta
    eta_plus[[index]] <- eta_plus[[index]] + step[[index]]
    eta_minus[[index]] <- eta_minus[[index]] - step[[index]]
    plus[[index]] <- build_covariance(eta_plus, design)$v
    minus[[index]] <- build_covariance(eta_minus, design)$v
    first[[index]] <- (plus[[index]] - minus[[index]]) / (2 * step[[index]])
    second[[index, index]] <- (
      plus[[index]] - 2 * base + minus[[index]]
    ) / step[[index]]^2
  }

  if (n_parameters > 1L) {
    for (row in 2:n_parameters) {
      for (column in seq_len(row - 1L)) {
        eta_pp <- eta_pm <- eta_mp <- eta_mm <- eta
        eta_pp[[row]] <- eta_pp[[row]] + step[[row]]
        eta_pp[[column]] <- eta_pp[[column]] + step[[column]]
        eta_pm[[row]] <- eta_pm[[row]] + step[[row]]
        eta_pm[[column]] <- eta_pm[[column]] - step[[column]]
        eta_mp[[row]] <- eta_mp[[row]] - step[[row]]
        eta_mp[[column]] <- eta_mp[[column]] + step[[column]]
        eta_mm[[row]] <- eta_mm[[row]] - step[[row]]
        eta_mm[[column]] <- eta_mm[[column]] - step[[column]]
        derivative <- (
          build_covariance(eta_pp, design)$v -
            build_covariance(eta_pm, design)$v -
            build_covariance(eta_mp, design)$v +
            build_covariance(eta_mm, design)$v
        ) / (4 * step[[row]] * step[[column]])
        second[[row, column]] <- derivative
        second[[column, row]] <- derivative
      }
    }
  }

  list(first = first, second = second)
}

kenward_roger_adjustment <- function(eta, eta_vcov, design, beta_vcov) {
  if (any(!is.finite(eta_vcov))) {
    cli::cli_abort(
      "Kenward-Roger inference requires an invertible likelihood Hessian."
    )
  }
  covariance <- build_covariance(eta, design)$v
  chol_v <- chol_factor(covariance)
  if (is.null(chol_v)) {
    cli::cli_abort(
      "Kenward-Roger inference requires a positive-definite covariance matrix."
    )
  }
  v_inverse <- chol_solve(chol_v, diag(nrow(covariance)))
  x <- design$x
  derivatives <- covariance_derivatives(eta, design)
  inverse_derivatives <- lapply(derivatives$first, function(derivative) {
    -v_inverse %*% derivative %*% v_inverse
  })
  p_matrices <- lapply(inverse_derivatives, function(derivative) {
    crossprod(x, derivative %*% x)
  })
  middle <- matrix(0, nrow(beta_vcov), ncol(beta_vcov))
  n_parameters <- length(eta)

  for (row in seq_len(n_parameters)) {
    for (column in seq_len(n_parameters)) {
      q_matrix <- crossprod(
        x,
        inverse_derivatives[[row]] %*%
          covariance %*%
          inverse_derivatives[[column]] %*%
          x
      )
      r_matrix <- crossprod(
        x,
        v_inverse %*%
          derivatives$second[[row, column]] %*%
          v_inverse %*%
          x
      )
      middle <- middle + eta_vcov[row, column] * (
        q_matrix -
          p_matrices[[row]] %*% beta_vcov %*% p_matrices[[column]] -
          0.25 * r_matrix
      )
    }
  }

  adjusted <- symmetrize(
    beta_vcov + 2 * beta_vcov %*% middle %*% beta_vcov
  )
  if (is.null(chol_factor(adjusted))) {
    message <- paste(
      "The Kenward-Roger adjusted fixed-effect covariance is not",
      "positive definite."
    )
    cli::cli_abort(message)
  }
  dimnames(adjusted) <- dimnames(beta_vcov)
  list(
    adjusted_vcov = adjusted,
    p_matrices = p_matrices,
    derivatives = derivatives
  )
}

kenward_roger_test_parameters <- function(object, contrast) {
  contrast <- if (is.null(dim(contrast))) {
    matrix(as.numeric(contrast), nrow = 1L)
  } else {
    as.matrix(contrast)
  }
  if (ncol(contrast) != length(object$coefficients)) {
    cli::cli_abort("The contrast has an incompatible number of columns.")
  }
  rank <- qr(contrast)$rank
  if (rank == 0L) {
    return(list(df = NA_real_, scale = NA_real_))
  }
  if (rank < nrow(contrast)) {
    row_qr <- qr(t(contrast))
    contrast <- contrast[row_qr$pivot[seq_len(rank)], , drop = FALSE]
  }

  v0 <- object$beta_vcov_model
  contrast_vcov <- symmetrize(contrast %*% v0 %*% t(contrast))
  inverse <- solve(contrast_vcov)
  m_matrix <- t(contrast) %*% inverse %*% contrast
  mv0 <- m_matrix %*% v0
  components <- lapply(object$kr$p_matrices, function(p_matrix) {
    mv0 %*% p_matrix %*% v0
  })
  a1 <- 0
  a2 <- 0
  w <- object$eta_vcov
  for (row in seq_along(components)) {
    for (column in seq_along(components)) {
      a1 <- a1 + w[row, column] *
        sum(diag(components[[row]])) *
        sum(diag(components[[column]]))
      a2 <- a2 + w[row, column] * sum(diag(
        components[[row]] %*% components[[column]]
      ))
    }
  }

  q <- nrow(contrast)
  if (!is.finite(a2) || abs(a2) < sqrt(.Machine$double.eps)) {
    return(list(df = Inf, scale = 1))
  }
  b <- (a1 + 6 * a2) / (2 * q)
  e_star <- 1 / (1 - a2 / q)
  g <- ((q + 1) * a1 - (q + 4) * a2) / ((q + 2) * a2)
  denominator <- 3 * q + 2 - 2 * g
  c1 <- g / denominator
  c2 <- (q - g) / denominator
  c3 <- (q + 2 - g) / denominator
  v_star <- 2 / q * (1 + c1 * b) /
    (1 - c2 * b)^2 / (1 - c3 * b)
  rho <- v_star / (2 * e_star^2)
  df <- 4 + (q + 2) / (q * rho - 1)
  scale <- df / (e_star * (df - 2))
  if (!is.finite(df) || df <= 2 || !is.finite(scale) || scale <= 0) {
    cli::cli_abort("Kenward-Roger degrees of freedom could not be computed.")
  }
  list(df = df, scale = scale)
}

beta_vcov_derivatives_at_eta <- function(eta, design, method = "Richardson") {
  p <- ncol(design$x)
  jacobian <- numDeriv::jacobian(
    func = function(parameters) {
      as.vector(beta_vcov_at_eta(parameters, design))
    },
    x = eta,
    method = method
  )
  derivatives <- lapply(seq_along(eta), function(index) {
    symmetrize(matrix(jacobian[, index], nrow = p, ncol = p))
  })
  names(derivatives) <- names(eta)
  derivatives
}

contrast_variance <- function(object, contrast, eta = object$eta) {
  beta_vcov <- beta_vcov_at_eta(eta, object$design)
  drop(contrast %*% beta_vcov %*% contrast)
}

satterthwaite_df <- function(object, contrast) {
  contrast <- as.numeric(contrast)
  variance <- drop(
    contrast %*% object$beta_vcov_model %*% contrast
  )
  invalid_variance <- !is.finite(variance) || variance <= 0
  if (invalid_variance || any(!is.finite(object$eta_vcov))) {
    return(NA_real_)
  }

  derivatives <- object$beta_vcov_derivatives %||%
    beta_vcov_derivatives_at_eta(
      object$eta,
      object$design,
      method = object$control$deriv_method
    )
  gradient <- vapply(derivatives, function(derivative) {
    drop(contrast %*% derivative %*% contrast)
  }, numeric(1L))
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
  if (object$ddf == "kenward-roger") {
    return(kenward_roger_test_parameters(object, contrast)$df)
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

  as_lmm_table(tibble::tibble(
    term = names(object$coefficients),
    estimate = unname(statistics[, "estimate"]),
    std.error = unname(statistics[, "std.error"]),
    statistic = unname(statistics[, "statistic"]),
    df = unname(statistics[, "df"]),
    p.value = unname(statistics[, "p.value"]),
    conf.low = unname(statistics[, "conf.low"]),
    conf.high = unname(statistics[, "conf.high"])
  ))
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
    cli::cli_warn(c(
      "A component denominator degree of freedom is at most 2.",
      "i" = paste(
        "The multi-degree-of-freedom test uses the conservative",
        "fallback of 2."
      )
    ))
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
  } else if (object$ddf == "kenward-roger") {
    parameters <- kenward_roger_test_parameters(
      object,
      projected_contrasts
    )
    denominator_df <- parameters$df
    statistic <- parameters$scale * statistic
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

  as_lmm_table(tibble::tibble(
    term = names(contrasts),
    num.df = unname(tests[, "num.df"]),
    den.df = unname(tests[, "den.df"]),
    statistic = unname(tests[, "statistic"]),
    p.value = unname(tests[, "p.value"])
  ))
}
