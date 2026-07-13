make_random_design <- function(random, data) {
  if (is.null(random)) {
    return(NULL)
  }

  parsed <- parse_bar_formula(random, "random")
  group <- make_group(data, parsed$group_vars, "random")
  random_frame <- stats::model.frame(
    parsed$terms_formula,
    data = data,
    na.action = stats::na.fail
  )
  z_base <- stats::model.matrix(parsed$terms_formula, random_frame)
  group_index <- as.integer(group)
  n_groups <- nlevels(group)
  n_terms <- ncol(z_base)
  n <- nrow(data)

  row_index <- rep(seq_len(n), each = n_terms)
  column_index <- as.vector(t(
    outer((group_index - 1L) * n_terms, seq_len(n_terms), "+")
  ))
  z <- Matrix::sparseMatrix(
    i = row_index,
    j = column_index,
    x = as.vector(t(z_base)),
    dims = c(n, n_groups * n_terms)
  )

  list(
    formula = random,
    parsed = parsed,
    group = group,
    group_index = group_index,
    group_levels = levels(group),
    group_label = paste(parsed$group_vars, collapse = ":"),
    z_base = z_base,
    z = z,
    term_names = colnames(z_base),
    n_groups = n_groups,
    n_terms = n_terms
  )
}

ordered_time <- function(x) {
  if (is.factor(x)) {
    levels_x <- levels(droplevels(x))
    return(list(index = match(as.character(x), levels_x), levels = levels_x))
  }
  if (inherits(x, c("Date", "POSIXct", "POSIXlt")) || is.numeric(x)) {
    levels_x <- sort(unique(x))
    return(list(index = match(x, levels_x), levels = as.character(levels_x)))
  }

  levels_x <- sort(unique(as.character(x)))
  list(index = match(as.character(x), levels_x), levels = levels_x)
}

make_repeated_design <- function(repeated, data, structure) {
  n <- nrow(data)
  if (is.null(repeated)) {
    if (structure != "id") {
      cli::cli_abort(
        "{.arg repeated} is required when {.arg structure} is not {.val id}."
      )
    }
    return(list(
      formula = NULL,
      group = factor(seq_len(n)),
      group_levels = as.character(seq_len(n)),
      group_label = "Observation",
      time_index = rep(1L, n),
      time_levels = "1",
      blocks = as.list(seq_len(n)),
      n_times = 1L
    ))
  }

  parsed <- parse_bar_formula(repeated, "repeated")
  group <- make_group(data, parsed$group_vars, "repeated")
  time_vars <- all.vars(parsed$terms_expr)
  if (length(time_vars) > 1L) {
    cli::cli_abort(
      "The left side of {.arg repeated} must contain one ordering variable."
    )
  }

  if (length(time_vars) == 0L) {
    time_index <- stats::ave(
      seq_len(n),
      group,
      FUN = function(index) seq_along(index)
    )
    time_levels <- as.character(seq_len(max(time_index)))
  } else {
    time_info <- ordered_time(data[[time_vars]])
    time_index <- time_info$index
    time_levels <- time_info$levels
  }

  key <- paste(as.integer(group), time_index, sep = "\r")
  if (anyDuplicated(key)) {
    cli::cli_abort(
      "Each repeated-measures group must have at most one observation per time."
    )
  }

  list(
    formula = repeated,
    parsed = parsed,
    group = group,
    group_levels = levels(group),
    group_label = paste(parsed$group_vars, collapse = ":"),
    time_variable = if (length(time_vars) == 1L) time_vars else NULL,
    time_index = as.integer(time_index),
    time_levels = time_levels,
    blocks = split(seq_len(n), group),
    n_times = length(time_levels)
  )
}

triangle_size <- function(n) {
  n * (n + 1L) / 2L
}

make_covariance_spec <- function(random_design, repeated_design, structure) {
  random_npar <- if (is.null(random_design)) {
    0L
  } else {
    triangle_size(random_design$n_terms)
  }

  q <- repeated_design$n_times
  residual_npar <- switch(structure,
    id = 1L,
    cs = 1L + as.integer(q > 1L),
    ar1 = 1L + as.integer(q > 1L),
    toep = q,
    un = triangle_size(q)
  )

  list(
    structure = structure,
    random_npar = random_npar,
    residual_npar = residual_npar,
    random_index = if (random_npar > 0L) seq_len(random_npar) else integer(),
    residual_index = random_npar + seq_len(residual_npar),
    npar = random_npar + residual_npar,
    n_times = q
  )
}

fill_cholesky <- function(eta, dimension) {
  out <- matrix(0, dimension, dimension)
  position <- 1L
  for (row in seq_len(dimension)) {
    for (column in seq_len(row)) {
      value <- eta[[position]]
      out[row, column] <- if (row == column) exp(value) else value
      position <- position + 1L
    }
  }
  out
}

pacf_to_acf <- function(partial_correlation) {
  order <- length(partial_correlation)
  if (order == 0L) {
    return(numeric())
  }

  ar <- partial_correlation[[1L]]
  if (order > 1L) {
    for (current_order in 2:order) {
      previous <- ar
      ar <- numeric(current_order)
      ar[[current_order]] <- partial_correlation[[current_order]]
      for (index in seq_len(current_order - 1L)) {
        ar[[index]] <- previous[[index]] -
          partial_correlation[[current_order]] *
            previous[[current_order - index]]
      }
    }
  }

  unname(stats::ARMAacf(ar = ar, lag.max = order)[-1L])
}

cs_correlation <- function(eta, n_times) {
  if (n_times <= 1L) {
    return(numeric())
  }
  epsilon <- 1e-7
  lower <- -1 / (n_times - 1L) + epsilon
  upper <- 1 - epsilon
  lower + (upper - lower) * stats::plogis(eta)
}

residual_covariance <- function(eta, spec) {
  structure <- spec$structure
  q <- spec$n_times

  if (structure == "un") {
    cholesky <- fill_cholesky(eta, q)
    return(tcrossprod(cholesky))
  }

  variance <- exp(eta[[1L]])
  correlation <- switch(structure,
    id = diag(q),
    cs = {
      rho <- if (q > 1L) cs_correlation(eta[[2L]], q) else 0
      matrix(rho, q, q) + diag(1 - rho, q)
    },
    ar1 = {
      rho <- if (q > 1L) tanh(eta[[2L]]) else 0
      rho^abs(outer(seq_len(q), seq_len(q), "-"))
    },
    toep = {
      partial <- if (q > 1L) tanh(eta[-1L]) else numeric()
      acf <- c(1, pacf_to_acf(partial))
      matrix(acf[abs(outer(seq_len(q), seq_len(q), "-")) + 1L], q, q)
    }
  )
  variance * correlation
}

random_covariance <- function(eta, random_design) {
  if (is.null(random_design)) {
    return(NULL)
  }
  cholesky <- fill_cholesky(eta, random_design$n_terms)
  tcrossprod(cholesky)
}

build_covariance <- function(eta, design) {
  spec <- design$covariance_spec
  random_design <- design$random
  repeated_design <- design$repeated
  n <- nrow(design$x)

  g <- NULL
  g_big <- NULL
  random_part <- matrix(0, n, n)
  if (!is.null(random_design)) {
    g <- random_covariance(eta[spec$random_index], random_design)
    g_big <- Matrix::bdiag(replicate(
      random_design$n_groups,
      g,
      simplify = FALSE
    ))
    z_matrix <- as.matrix(random_design$z)
    random_part <- as.matrix(z_matrix %*% g_big %*% t(z_matrix))
  }

  residual_base <- residual_covariance(eta[spec$residual_index], spec)
  if (spec$structure == "id") {
    residual <- diag(residual_base[[1L]], n)
  } else {
    residual <- matrix(0, n, n)
    for (rows in repeated_design$blocks) {
      times <- repeated_design$time_index[rows]
      residual[rows, rows] <- residual_base[times, times, drop = FALSE]
    }
  }

  list(
    v = symmetrize(random_part + residual),
    g = g,
    g_big = g_big,
    r = residual,
    residual_base = residual_base
  )
}

initial_eta <- function(y, design) {
  spec <- design$covariance_spec
  response_variance <- stats::var(y)
  if (!is.finite(response_variance) || response_variance <= 0) {
    response_variance <- 1
  }

  out <- numeric(spec$npar)
  if (spec$random_npar > 0L) {
    k <- design$random$n_terms
    position <- 1L
    for (row in seq_len(k)) {
      for (column in seq_len(row)) {
        out[[position]] <- if (row == column) {
          log(sqrt(0.25 * response_variance / k))
        } else {
          0
        }
        position <- position + 1L
      }
    }
  }

  residual_index <- spec$residual_index
  if (spec$structure == "un") {
    q <- spec$n_times
    position <- 1L
    for (row in seq_len(q)) {
      for (column in seq_len(row)) {
        out[residual_index[[position]]] <- if (row == column) {
          log(sqrt(0.75 * response_variance))
        } else {
          0
        }
        position <- position + 1L
      }
    }
  } else {
    out[residual_index[[1L]]] <- log(0.75 * response_variance)
    if (length(residual_index) > 1L) {
      q <- spec$n_times
      if (spec$structure == "cs") {
        lower <- -1 / (q - 1L) + 1e-7
        upper <- 1 - 1e-7
        probability <- (0.1 - lower) / (upper - lower)
        out[residual_index[[2L]]] <- stats::qlogis(probability)
      } else if (spec$structure == "ar1") {
        out[residual_index[[2L]]] <- atanh(0.1)
      } else {
        out[residual_index[-1L]] <- c(atanh(0.1), rep(0, q - 2L))
      }
    }
  }
  out
}

covariance_natural <- function(eta, design) {
  spec <- design$covariance_spec
  values <- numeric()

  if (!is.null(design$random)) {
    g <- random_covariance(eta[spec$random_index], design$random)
    terms <- design$random$term_names
    random_values <- stats::setNames(diag(g), paste0("random.var.", terms))
    if (length(terms) > 1L) {
      correlations <- stats::cov2cor(g)
      for (row in 2:length(terms)) {
        for (column in seq_len(row - 1L)) {
          name <- paste0("random.cor.", terms[[row]], ".", terms[[column]])
          random_values[[name]] <- correlations[row, column]
        }
      }
    }
    values <- c(values, random_values)
  }

  residual <- residual_covariance(eta[spec$residual_index], spec)
  time_names <- design$repeated$time_levels
  if (spec$structure == "un") {
    residual_values <- stats::setNames(
      diag(residual),
      paste0("residual.var.", time_names)
    )
    correlations <- stats::cov2cor(residual)
    if (length(time_names) > 1L) {
      for (row in 2:length(time_names)) {
        for (column in seq_len(row - 1L)) {
          name <- paste0(
            "residual.cor.", time_names[[row]], ".", time_names[[column]]
          )
          residual_values[[name]] <- correlations[row, column]
        }
      }
    }
  } else {
    residual_values <- c("residual.var" = residual[[1L]])
    if (spec$structure == "cs" && spec$n_times > 1L) {
      residual_values[["residual.cor"]] <- residual[1L, 2L] / residual[[1L]]
    }
    if (spec$structure == "ar1" && spec$n_times > 1L) {
      residual_values[["residual.cor"]] <-
        residual[1L, 2L] / residual[[1L]]
    }
    if (spec$structure == "toep" && spec$n_times > 1L) {
      for (lag in seq_len(spec$n_times - 1L)) {
        residual_values[[paste0("residual.cor.lag", lag)]] <-
          residual[1L, lag + 1L] / residual[[1L]]
      }
    }
  }

  c(values, residual_values)
}
