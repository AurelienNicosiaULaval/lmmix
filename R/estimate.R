run_covariance_optimizer <- function(
  start,
  objective,
  optimizer,
  optim_method,
  lower,
  upper,
  control
) {
  if (optimizer == "nlminb") {
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
    list(
      eta = optimization$par,
      code = optimization$convergence,
      message = optimization$message,
      iterations = optimization$iterations,
      evaluations = optimization$evaluations,
      optimizer = "nlminb",
      method = NA_character_
    )
  } else {
    optimization <- stats::optim(
      par = start,
      fn = objective,
      method = optim_method,
      lower = if (optim_method == "L-BFGS-B") lower else -Inf,
      upper = if (optim_method == "L-BFGS-B") upper else Inf,
      control = list(
        maxit = control$max_iter,
        reltol = control$rel_tol
      )
    )
    list(
      eta = optimization$par,
      code = optimization$convergence,
      message = optimization$message %||% "",
      iterations = unname(optimization$counts[["function"]]),
      evaluations = optimization$counts,
      optimizer = "optim",
      method = optim_method
    )
  }
}

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

  strategies <- if (control$optimizer == "auto") {
    list(
      c(optimizer = "nlminb", method = NA_character_),
      c(optimizer = "optim", method = "BFGS"),
      c(optimizer = "optim", method = "L-BFGS-B")
    )
  } else {
    list(c(
      optimizer = control$optimizer,
      method = if (control$optimizer == "optim") {
        control$optim_method
      } else {
        NA_character_
      }
    ))
  }

  max_attempts <- 1L + control$max_restarts
  attempts <- vector("list", max_attempts)
  for (attempt in seq_len(max_attempts)) {
    strategy <- strategies[[1L + (attempt - 1L) %% length(strategies)]]
    attempt_start <- start
    if (attempt > 1L && control$restart_scale > 0) {
      perturbation <- control$restart_scale * sin(
        seq_len(npar) * (attempt - 1L)
      )
      attempt_start <- pmin(upper, pmax(lower, start + perturbation))
    }
    result <- tryCatch(
      run_covariance_optimizer(
        attempt_start,
        objective,
        optimizer = strategy[["optimizer"]],
        optim_method = strategy[["method"]],
        lower = lower,
        upper = upper,
        control = control
      ),
      error = function(cnd) {
        list(
          eta = attempt_start,
          code = 999L,
          message = conditionMessage(cnd),
          iterations = NA_integer_,
          evaluations = NA,
          optimizer = strategy[["optimizer"]],
          method = strategy[["method"]]
        )
      }
    )
    result$objective <- objective(result$eta)
    result$attempt <- attempt
    attempts[[attempt]] <- result
    if (result$code == 0L && is.finite(result$objective)) {
      attempts <- attempts[seq_len(attempt)]
      break
    }
  }

  converged <- vapply(attempts, function(x) x$code == 0L, logical(1L))
  candidates <- if (any(converged)) which(converged) else seq_along(attempts)
  selected_index <- candidates[[which.min(vapply(
    attempts[candidates],
    `[[`,
    numeric(1L),
    "objective"
  ))]]
  selected <- attempts[[selected_index]]
  eta <- selected$eta

  hessian <- numDeriv::hessian(objective, eta, method = control$deriv_method)
  hessian_result <- invert_hessian(hessian)
  if (selected$code != 0L) {
    cli::cli_warn(c(
      "The optimizer did not report convergence.",
      "i" = "Code {selected$code}: {selected$message}"
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
    code = selected$code,
    message = selected$message,
    iterations = selected$iterations,
    evaluations = selected$evaluations,
    optimizer = selected$optimizer,
    attempts = tibble::tibble(
      attempt = vapply(attempts, `[[`, integer(1L), "attempt"),
      optimizer = vapply(attempts, `[[`, character(1L), "optimizer"),
      method = vapply(attempts, `[[`, character(1L), "method"),
      convergence = vapply(attempts, `[[`, integer(1L), "code"),
      objective = vapply(attempts, `[[`, numeric(1L), "objective"),
      selected = seq_along(attempts) == selected_index
    )
  )
}

estimate_blup <- function(gls, design) {
  random_design <- design$random
  if (is.null(random_design)) {
    return(NULL)
  }

  effects <- lapply(seq_along(random_design$terms), function(index) {
    term <- random_design$terms[[index]]
    label <- names(random_design$terms)[[index]]
    u <- drop(
      gls$covariance$g_big[[label]] %*%
        Matrix::t(term$z) %*%
        gls$v_inv_residual
    )
    matrix(
      u,
      nrow = term$n_groups,
      ncol = term$n_terms,
      byrow = TRUE,
      dimnames = list(term$group_levels, term$term_names)
    )
  })
  names(effects) <- names(random_design$terms)
  effects
}

covariance_metadata <- function(names, design) {
  rows <- lapply(names, function(name) {
    pieces <- strsplit(name, ".", fixed = TRUE)[[1L]]
    if (pieces[[1L]] == "random") {
      component <- pieces[[2L]]
      multiple <- length(design$random$terms) > 1L
      random_index <- if (multiple) as.integer(pieces[[3L]]) else 1L
      term_start <- if (multiple) 4L else 3L
      random_group <- design$random$terms[[random_index]]$group_label
      separator <- if (component == "cor") ", " else "."
      term <- paste(pieces[term_start:length(pieces)], collapse = separator)
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
  as_lmm_table(tibble::add_column(
    metadata,
    estimate = unname(natural),
    std.error = standard_error,
    .after = "component"
  ))
}
