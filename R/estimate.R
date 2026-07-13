fit_covariance_parameters <- function(design, method, control) {
  objective <- function(eta) {
    negative_log_likelihood(eta, design = design, method = method)
  }
  npar <- design$covariance_spec$npar
  start <- control$initial %||% initial_eta(design$y, design)
  if (!is.numeric(start) || length(start) != npar || any(!is.finite(start))) {
    cli::cli_abort(
      "{.arg control$initial} must contain {npar} finite numeric value{?s}."
    )
  }

  lower <- rep(control$lower, length.out = npar)
  upper <- rep(control$upper, length.out = npar)
  if (any(start < lower | start > upper)) {
    cli::cli_abort(
      "{.arg control$initial} lies outside the optimization bounds."
    )
  }

  if (control$optimizer == "nlminb") {
    optimization <- stats::nlminb(
      start = start,
      objective = objective,
      lower = lower,
      upper = upper,
      control = list(
        iter.max = control$max_iter,
        eval.max = max(200L, 2L * control$max_iter),
        rel.tol = control$rel_tol,
        x.tol = control$x_tol
      )
    )
    eta <- optimization$par
    convergence_code <- optimization$convergence
    convergence_message <- optimization$message
    iterations <- optimization$iterations
    evaluations <- optimization$evaluations
  } else {
    optimization <- stats::optim(
      par = start,
      fn = objective,
      method = control$optim_method,
      lower = if (control$optim_method == "L-BFGS-B") lower else -Inf,
      upper = if (control$optim_method == "L-BFGS-B") upper else Inf,
      control = list(
        maxit = control$max_iter,
        reltol = control$rel_tol
      )
    )
    eta <- optimization$par
    convergence_code <- optimization$convergence
    convergence_message <- optimization$message %||% ""
    iterations <- unname(optimization$counts[["function"]])
    evaluations <- optimization$counts
  }

  hessian <- numDeriv::hessian(objective, eta, method = control$deriv_method)
  hessian_result <- invert_hessian(hessian)
  if (convergence_code != 0L) {
    cli::cli_warn(c(
      "The optimizer did not report convergence.",
      "i" = "Code {convergence_code}: {convergence_message}"
    ))
  }
  if (!hessian_result$positive_definite) {
    message <- paste(
      "The likelihood Hessian is not positive definite;",
      "some inference may be unavailable."
    )
    cli::cli_warn(message)
  }

  list(
    eta = eta,
    objective = objective(eta),
    hessian = hessian,
    eta_vcov = hessian_result$vcov,
    hessian_positive_definite = hessian_result$positive_definite,
    hessian_eigenvalues = hessian_result$eigenvalues,
    code = convergence_code,
    message = convergence_message,
    iterations = iterations,
    evaluations = evaluations,
    optimizer = control$optimizer
  )
}

estimate_blup <- function(gls, design) {
  random_design <- design$random
  if (is.null(random_design)) {
    return(NULL)
  }

  u <- drop(
    gls$covariance$g_big %*%
      t(as.matrix(random_design$z)) %*%
      gls$v_inv_residual
  )
  matrix(
    u,
    nrow = random_design$n_groups,
    ncol = random_design$n_terms,
    byrow = TRUE,
    dimnames = list(
      random_design$group_levels,
      random_design$term_names
    )
  )
}

covariance_metadata <- function(names, design) {
  random_group <- if (is.null(design$random)) {
    NA_character_
  } else {
    design$random$group_label
  }

  rows <- lapply(names, function(name) {
    pieces <- strsplit(name, ".", fixed = TRUE)[[1L]]
    if (pieces[[1L]] == "random") {
      component <- pieces[[2L]]
      separator <- if (component == "cor") ", " else "."
      term <- paste(pieces[-c(1L, 2L)], collapse = separator)
      list(group = random_group, term = term, component = component)
    } else {
      component <- pieces[[2L]]
      term <- if (length(pieces) > 2L) {
        paste(pieces[-c(1L, 2L)], collapse = ".")
      } else {
        design$covariance_spec$structure
      }
      list(group = "Residual", term = term, component = component)
    }
  })

  tibble::tibble(
    group = vapply(rows, `[[`, character(1L), "group"),
    term = vapply(rows, `[[`, character(1L), "term"),
    component = vapply(rows, `[[`, character(1L), "component")
  )
}

estimate_covariance_components <- function(eta, eta_vcov, design) {
  natural <- covariance_natural(eta, design)
  standard_error <- rep(NA_real_, length(natural))
  if (all(is.finite(eta_vcov))) {
    jacobian <- numDeriv::jacobian(
      func = covariance_natural,
      x = eta,
      design = design
    )
    natural_vcov <- symmetrize(jacobian %*% eta_vcov %*% t(jacobian))
    standard_error <- sqrt(pmax(diag(natural_vcov), 0))
  }

  metadata <- covariance_metadata(names(natural), design)
  dot_names(tibble::add_column(
    metadata,
    estimate = unname(natural),
    std.error = standard_error,
    .after = "component"
  ))
}
