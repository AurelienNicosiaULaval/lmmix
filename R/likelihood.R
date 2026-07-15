gls_at_eta_dense <- function(eta, design) {
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

gls_at_eta_blockwise <- function(eta, design, keep_covariance = FALSE) {
  x <- design$x
  y <- design$y
  p <- ncol(x)
  parts <- covariance_parts(eta, design)
  information <- matrix(0, p, p)
  score <- numeric(p)
  block_factors <- vector("list", length(design$likelihood_blocks))
  logdet_v <- 0

  for (index in seq_along(design$likelihood_blocks)) {
    rows <- design$likelihood_blocks[[index]]
    block <- build_covariance_block(eta, design, rows, parts = parts)
    chol_v <- chol_factor(block$v)
    if (is.null(chol_v)) {
      return(NULL)
    }
    x_block <- x[rows, , drop = FALSE]
    y_block <- y[rows]
    v_inv_x <- chol_solve(chol_v, x_block)
    v_inv_y <- drop(chol_solve(chol_v, y_block))
    information <- information + crossprod(x_block, v_inv_x)
    score <- score + drop(crossprod(x_block, v_inv_y))
    logdet_v <- logdet_v + chol_logdet(chol_v)
    block_factors[[index]] <- list(rows = rows, chol = chol_v)
  }

  chol_information <- chol_factor(information)
  if (is.null(chol_information)) {
    return(NULL)
  }
  beta <- drop(chol_solve(chol_information, score))
  names(beta) <- colnames(x)
  beta_vcov <- chol_solve(chol_information, diag(p))
  dimnames(beta_vcov) <- list(colnames(x), colnames(x))
  residual <- y - drop(x %*% beta)
  v_inv_residual <- numeric(length(y))
  quadratic <- 0
  for (block in block_factors) {
    solved <- drop(chol_solve(block$chol, residual[block$rows]))
    v_inv_residual[block$rows] <- solved
    quadratic <- quadratic + sum(residual[block$rows] * solved)
  }

  list(
    beta = beta,
    beta_vcov = symmetrize(beta_vcov),
    residual = residual,
    v_inv_residual = v_inv_residual,
    quadratic = quadratic,
    logdet_v = logdet_v,
    logdet_information = chol_logdet(chol_information),
    chol_v = lapply(block_factors, `[[`, "chol"),
    covariance = if (isTRUE(keep_covariance)) {
      build_covariance(eta, design)
    } else {
      NULL
    },
    blockwise = TRUE
  )
}

gls_at_eta <- function(eta, design, keep_covariance = FALSE) {
  blocks <- design$likelihood_blocks %||% list(seq_len(nrow(design$x)))
  if (length(blocks) <= 1L) {
    return(gls_at_eta_dense(eta, design))
  }
  gls_at_eta_blockwise(eta, design, keep_covariance = keep_covariance)
}

nll_from_gls <- function(fit, n, p, method) {
  if (is.null(fit)) {
    return(.Machine$double.xmax / 1000)
  }

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
  nll_from_gls(fit, n, p, method)
}

beta_vcov_at_eta <- function(eta, design) {
  result <- gls_at_eta(eta, design)
  if (is.null(result)) {
    cli::cli_abort("The covariance matrix is not positive definite.")
  }
  result$beta_vcov
}
