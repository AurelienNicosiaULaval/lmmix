make_random_design_term <- function(random, data) {
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

make_random_design <- function(random, data) {
  formulas <- normalize_random_formulas(random)
  if (is.null(formulas)) {
    return(NULL)
  }

  terms <- lapply(formulas, make_random_design_term, data = data)
  supplied_names <- names(formulas)
  labels <- vapply(terms, `[[`, character(1L), "group_label")
  if (!is.null(supplied_names)) {
    use_name <- nzchar(supplied_names)
    labels[use_name] <- supplied_names[use_name]
  }
  labels <- make.unique(labels)
  names(terms) <- labels
  names(formulas) <- labels

  list(
    formulas = formulas,
    terms = terms,
    labels = labels,
    z = do.call(cbind, lapply(terms, `[[`, "z"))
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

make_likelihood_blocks <- function(random_design, repeated_design, n) {
  parent <- seq_len(n)

  find_root <- function(index) {
    root <- index
    while (parent[[root]] != root) {
      root <- parent[[root]]
    }
    while (parent[[index]] != index) {
      next_index <- parent[[index]]
      parent[[index]] <<- root
      index <- next_index
    }
    root
  }

  merge_rows <- function(rows) {
    if (length(rows) <= 1L) {
      return(invisible(NULL))
    }
    root <- find_root(rows[[1L]])
    for (index in rows[-1L]) {
      other <- find_root(index)
      if (other != root) {
        parent[[other]] <<- root
      }
    }
    invisible(NULL)
  }

  lapply(repeated_design$blocks, merge_rows)
  if (!is.null(random_design)) {
    for (term in random_design$terms) {
      lapply(split(seq_len(n), term$group_index), merge_rows)
    }
  }

  roots <- vapply(seq_len(n), find_root, integer(1L))
  unname(split(seq_len(n), match(roots, unique(roots))))
}

triangle_size <- function(n) {
  n * (n + 1L) / 2L
}

parse_covariance_structure <- function(structure) {
  if (!is.character(structure) || length(structure) != 1L || is.na(structure)) {
    cli::cli_abort("{.arg structure} must be a single character value.")
  }
  value <- tolower(gsub("[[:space:]]+", "", structure))
  if (value %in% c("id", "cs", "ar1", "toep", "un")) {
    return(list(name = value, order = NULL, label = value))
  }

  match <- regexec("^toep\\(([1-9][0-9]*)\\)$", value)
  pieces <- regmatches(value, match)[[1L]]
  if (length(pieces) == 2L) {
    order <- as.integer(pieces[[2L]])
    return(list(
      name = "toep",
      order = order,
      label = paste0("toep(", order, ")")
    ))
  }

  cli::cli_abort(
    paste(
      "{.arg structure} must be one of {.val id}, {.val cs}, {.val ar1},",
      "{.val toep}, {.code toep(k)}, or {.val un}."
    )
  )
}

make_covariance_spec <- function(
  random_design,
  repeated_design,
  structure,
  covariance_order = NULL
) {
  random_npar_by_term <- if (is.null(random_design)) {
    integer()
  } else {
    vapply(
      random_design$terms,
      function(term) triangle_size(term$n_terms),
      numeric(1L)
    )
  }
  random_npar <- sum(random_npar_by_term)
  random_index <- if (random_npar == 0L) {
    list()
  } else {
    ends <- cumsum(random_npar_by_term)
    starts <- ends - random_npar_by_term + 1L
    Map(seq.int, starts, ends)
  }
  names(random_index) <- names(random_npar_by_term)

  q <- repeated_design$n_times
  explicit_order <- !is.null(covariance_order)
  if (structure == "toep") {
    covariance_order <- covariance_order %||% q
    if (covariance_order > q) {
      cli::cli_abort(
        "The order in {.code toep(k)} cannot exceed the {q} time levels."
      )
    }
  }
  residual_npar <- switch(structure,
    id = 1L,
    cs = 1L + as.integer(q > 1L),
    ar1 = 1L + as.integer(q > 1L),
    toep = covariance_order,
    un = triangle_size(q)
  )

  list(
    structure = structure,
    random_npar = random_npar,
    random_npar_by_term = random_npar_by_term,
    residual_npar = residual_npar,
    random_index = random_index,
    residual_index = random_npar + seq_len(residual_npar),
    npar = random_npar + residual_npar,
    n_times = q,
    covariance_order = if (structure == "toep") {
      as.integer(covariance_order)
    } else {
      NA_integer_
    },
    structure_label = if (structure == "toep" && explicit_order) {
      paste0("toep(", covariance_order, ")")
    } else {
      structure
    }
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
      order <- spec$covariance_order
      partial <- if (order > 1L) tanh(eta[-1L]) else numeric()
      acf <- c(
        1,
        pacf_to_acf(partial),
        rep(0, q - order)
      )
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

covariance_parts <- function(eta, design) {
  spec <- design$covariance_spec
  random_design <- design$random

  g <- list()
  g_big <- list()
  if (!is.null(random_design)) {
    for (index in seq_along(random_design$terms)) {
      term <- random_design$terms[[index]]
      label <- names(random_design$terms)[[index]]
      term_g <- random_covariance(eta[spec$random_index[[index]]], term)
      term_g_big <- Matrix::bdiag(replicate(
        term$n_groups,
        term_g,
        simplify = FALSE
      ))
      g[[label]] <- term_g
      g_big[[label]] <- term_g_big
    }
  }

  residual_base <- residual_covariance(eta[spec$residual_index], spec)

  list(
    g = if (length(g) == 0L) NULL else g,
    g_big = if (length(g_big) == 0L) NULL else g_big,
    residual_base = residual_base
  )
}

build_covariance_block <- function(eta, design, rows, parts = NULL) {
  parts <- parts %||% covariance_parts(eta, design)
  repeated_design <- design$repeated
  n <- length(rows)

  random_part <- matrix(0, n, n)
  if (!is.null(design$random)) {
    for (label in names(design$random$terms)) {
      z <- design$random$terms[[label]]$z[rows, , drop = FALSE]
      random_part <- random_part + as.matrix(
        z %*% parts$g_big[[label]] %*% Matrix::t(z)
      )
    }
  }

  if (design$covariance_spec$structure == "id") {
    residual <- diag(parts$residual_base[[1L]], n)
  } else {
    residual <- matrix(0, n, n)
    local_blocks <- split(seq_along(rows), repeated_design$group[rows])
    for (positions in local_blocks) {
      global_rows <- rows[positions]
      times <- repeated_design$time_index[global_rows]
      residual[positions, positions] <- parts$residual_base[
        times,
        times,
        drop = FALSE
      ]
    }
  }
  residual_scale <- 1 / sqrt(design$weights[rows])
  residual <- residual * tcrossprod(residual_scale)

  list(
    v = symmetrize(random_part + residual),
    r = residual
  )
}

build_covariance <- function(eta, design) {
  parts <- covariance_parts(eta, design)
  block <- build_covariance_block(
    eta,
    design,
    rows = seq_len(nrow(design$x)),
    parts = parts
  )

  list(
    v = block$v,
    g = parts$g,
    g_big = parts$g_big,
    r = block$r,
    residual_base = parts$residual_base
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
    n_random <- length(design$random$terms)
    for (index in seq_along(design$random$terms)) {
      k <- design$random$terms[[index]]$n_terms
      term_index <- spec$random_index[[index]]
      position <- 1L
      for (row in seq_len(k)) {
        for (column in seq_len(row)) {
          out[term_index[[position]]] <- if (row == column) {
            log(sqrt(0.25 * response_variance / (n_random * k)))
          } else {
            0
          }
          position <- position + 1L
        }
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
        correlation_count <- length(residual_index) - 1L
        out[residual_index[-1L]] <- c(
          atanh(0.1),
          rep(0, correlation_count - 1L)
        )
      }
    }
  }
  out
}

covariance_natural <- function(eta, design) {
  spec <- design$covariance_spec
  values <- numeric()

  if (!is.null(design$random)) {
    multiple <- length(design$random$terms) > 1L
    for (index in seq_along(design$random$terms)) {
      random_term <- design$random$terms[[index]]
      g <- random_covariance(eta[spec$random_index[[index]]], random_term)
      terms <- random_term$term_names
      prefix <- if (multiple) {
        paste0("random.var.", index, ".")
      } else {
        "random.var."
      }
      random_values <- stats::setNames(diag(g), paste0(prefix, terms))
      if (length(terms) > 1L) {
        correlations <- stats::cov2cor(g)
        for (row in 2:length(terms)) {
          for (column in seq_len(row - 1L)) {
            prefix <- if (multiple) {
              paste0("random.cor.", index, ".")
            } else {
              "random.cor."
            }
            name <- paste0(prefix, terms[[row]], ".", terms[[column]])
            random_values[[name]] <- correlations[row, column]
          }
        }
      }
      values <- c(values, random_values)
    }
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
      for (lag in seq_len(spec$covariance_order - 1L)) {
        residual_values[[paste0("residual.cor.lag", lag)]] <-
          residual[1L, lag + 1L] / residual[[1L]]
      }
    }
  }

  c(values, residual_values)
}
