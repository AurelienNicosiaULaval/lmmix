gls_at_eta <- function(eta, design) {
  covariance <- build_covariance(eta, design)
  chol_v <- chol_factor(covariance$v)
  if (is.null(chol_v)) {
    return(NULL)
  }

  x <- design$x
  y <- design$y
  v_inv_x <- chol_solve(chol_v, x)
  v_inv_y <- drop(chol_solve(chol_v, y))
  information <- crossprod(x, v_inv_x)
  chol_information <- chol_factor(information)
  if (is.null(chol_information)) {
    return(NULL)
  }

  beta <- drop(chol_solve(chol_information, crossprod(x, v_inv_y)))
  names(beta) <- colnames(x)
  beta_vcov <- chol_solve(chol_information, diag(ncol(x)))
  dimnames(beta_vcov) <- list(colnames(x), colnames(x))
  residual <- y - drop(x %*% beta)
  v_inv_residual <- drop(chol_solve(chol_v, residual))

  list(
    beta = beta,
    beta_vcov = symmetrize(beta_vcov),
    residual = residual,
    v_inv_residual = v_inv_residual,
    quadratic = sum(residual * v_inv_residual),
    logdet_v = chol_logdet(chol_v),
    logdet_information = chol_logdet(chol_information),
    chol_v = chol_v,
    covariance = covariance
  )
}

negative_log_likelihood <- function(eta, design, method) {
  fit <- tryCatch(
    gls_at_eta(eta, design),
    error = function(cnd) NULL
  )
  if (is.null(fit)) {
    return(.Machine$double.xmax / 1000)
  }

  n <- nrow(design$x)
  p <- ncol(design$x)
  if (method == "REML") {
    0.5 * (
      (n - p) * log(2 * pi) +
        fit$logdet_v +
        fit$logdet_information +
        fit$quadratic
    )
  } else {
    0.5 * (
      n * log(2 * pi) +
        fit$logdet_v +
        fit$quadratic
    )
  }
}

beta_vcov_at_eta <- function(eta, design) {
  result <- gls_at_eta(eta, design)
  if (is.null(result)) {
    cli::cli_abort("The covariance matrix is not positive definite.")
  }
  result$beta_vcov
}
