parse_lsmeans_specs <- function(specs) {
  pairwise <- FALSE
  if (is.character(specs)) {
    variables <- specs
  } else if (inherits(specs, "formula")) {
    if (length(specs) == 3L) {
      left <- paste(deparse(specs[[2L]]), collapse = "")
      pairwise <- identical(tolower(left), "pairwise")
      if (!pairwise) {
        cli::cli_abort(
          "The left side of {.arg specs} must be {.code pairwise}."
        )
      }
      variables <- all.vars(specs[[3L]])
    } else {
      variables <- all.vars(specs[[2L]])
    }
  } else {
    cli::cli_abort("{.arg specs} must be a character vector or formula.")
  }

  variables <- unique(variables)
  if (length(variables) == 0L) {
    cli::cli_abort("{.arg specs} must identify at least one variable.")
  }
  list(variables = variables, pairwise = pairwise)
}

reference_values <- function(x, requested, at = NULL) {
  if (!is.null(at)) {
    return(at)
  }
  if (is.factor(x)) {
    return(levels(x))
  }
  if (is.character(x)) {
    return(sort(unique(x)))
  }
  if (is.logical(x)) {
    return(c(FALSE, TRUE))
  }
  if (requested) {
    return(sort(unique(x)))
  }
  if (is.numeric(x)) {
    return(mean(x))
  }

  cli::cli_abort(
    "Cannot construct a reference value for class {.cls {class(x)}}."
  )
}

make_reference_grid <- function(object, specs, at) {
  variables <- all.vars(stats::delete.response(object$terms))
  missing_vars <- setdiff(c(variables, specs, names(at)), names(object$data))
  if (length(missing_vars) > 0L) {
    cli::cli_abort(
      "Unknown reference-grid variable{?s}: {.field {missing_vars}}."
    )
  }
  if (length(setdiff(specs, variables)) > 0L) {
    cli::cli_abort(
      "Every variable in {.arg specs} must occur in the fixed-effects formula."
    )
  }

  values <- lapply(variables, function(variable) {
    reference_values(
      object$data[[variable]],
      requested = variable %in% specs,
      at = at[[variable]]
    )
  })
  names(values) <- variables
  grid <- do.call(
    expand.grid,
    c(values, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  )

  for (variable in variables) {
    original <- object$data[[variable]]
    if (is.factor(original)) {
      grid[[variable]] <- factor(
        grid[[variable]],
        levels = levels(original),
        ordered = is.ordered(original)
      )
    }
  }
  grid
}

lsmeans_contrast_matrix <- function(object, specs, at = list()) {
  if (!is.list(at) || (length(at) > 0L && is.null(names(at)))) {
    cli::cli_abort("{.arg at} must be a named list.")
  }
  grid <- make_reference_grid(object, specs, at)
  x_grid <- model_matrix_newdata(object, grid)
  coefficient_names <- names(object$coefficients)
  if (!setequal(colnames(x_grid), coefficient_names)) {
    cli::cli_abort("The reference-grid model matrix is not estimable.")
  }
  x_grid <- x_grid[, coefficient_names, drop = FALSE]

  combinations <- unique(grid[specs])
  contrast <- matrix(
    0,
    nrow = nrow(combinations),
    ncol = length(coefficient_names),
    dimnames = list(NULL, coefficient_names)
  )
  for (index in seq_len(nrow(combinations))) {
    selected <- rep(TRUE, nrow(grid))
    for (variable in specs) {
      selected <- selected & grid[[variable]] == combinations[[variable]][index]
    }
    contrast[index, ] <- colMeans(x_grid[selected, , drop = FALSE])
  }

  list(
    combinations = combinations,
    contrast = contrast,
    grid = grid
  )
}

format_lsmean_label <- function(row, specs) {
  values <- vapply(
    specs,
    function(variable) as.character(row[[variable]]),
    character(1L)
  )
  if (length(specs) == 1L) {
    values
  } else {
    paste(paste0(specs, "=", values), collapse = ", ")
  }
}

lsmeans_estimates <- function(object, information, level) {
  statistics <- lapply(seq_len(nrow(information$contrast)), function(index) {
    contrast_statistics(object, information$contrast[index, ], level = level)
  })
  statistics <- do.call(rbind, statistics)

  out <- tibble::as_tibble(information$combinations)
  out$estimate <- statistics[, "estimate"]
  out$std.error <- statistics[, "std.error"]
  out$df <- statistics[, "df"]
  out$statistic <- statistics[, "statistic"]
  out$p.value <- statistics[, "p.value"]
  out$conf.low <- statistics[, "conf.low"]
  out$conf.high <- statistics[, "conf.high"]
  class(out) <- c("lmm_lsmeans", class(out))
  attr(out, "contrast_matrix") <- information$contrast
  out
}

lsmeans_pairs <- function(object, information, specs, level, adjust) {
  n_means <- nrow(information$contrast)
  if (n_means < 2L) {
    cli::cli_abort("At least two LS-means are required for pairwise contrasts.")
  }
  pairs <- utils::combn(seq_len(n_means), 2L)
  statistics <- lapply(seq_len(ncol(pairs)), function(index) {
    first <- pairs[1L, index]
    second <- pairs[2L, index]
    contrast <- information$contrast[first, ] - information$contrast[second, ]
    contrast_statistics(object, contrast, level = level)
  })
  statistics <- do.call(rbind, statistics)

  labels <- vapply(seq_len(nrow(information$combinations)), function(index) {
    format_lsmean_label(information$combinations[index, , drop = FALSE], specs)
  }, character(1L))
  out <- tibble::tibble(
    contrast = paste0(labels[pairs[1L, ]], " - ", labels[pairs[2L, ]]),
    estimate = statistics[, "estimate"],
    std.error = statistics[, "std.error"],
    df = statistics[, "df"],
    statistic = statistics[, "statistic"],
    p.value = stats::p.adjust(statistics[, "p.value"], method = adjust),
    conf.low = statistics[, "conf.low"],
    conf.high = statistics[, "conf.high"]
  )
  class(out) <- c("lmm_lsmeans_contrasts", class(out))
  out
}

#' Estimated marginal means for an `lmm` model
#'
#' Marginal means use equal weights over nuisance-factor levels. Numeric
#' covariates not listed in `specs` are held at their observed means.
#' Multiplicity adjustment applies to pairwise p-values; confidence intervals
#' remain pointwise.
#'
#' @param object An `lmm` object.
#' @param specs Variables defining the marginal means, supplied as a character
#'   vector or one-sided formula. Use `pairwise ~ factor` to request
#'   differences.
#' @param pairwise Whether to return pairwise differences as well as means.
#' @param at Named list overriding reference-grid values.
#' @param level Confidence level.
#' @param adjust Multiplicity adjustment passed to [stats::p.adjust()].
#' @param ... Reserved for future extensions.
#'
#' @return A tibble, or a list with `lsmeans` and `contrasts` tibbles when
#'   pairwise comparisons are requested.
#' @export
lsmeans <- function(
  object,
  specs,
  pairwise = FALSE,
  at = list(),
  level = 0.95,
  adjust = "none",
  ...
) {
  if (!inherits(object, "lmm")) {
    cli::cli_abort("{.arg object} must inherit from {.cls lmm}.")
  }
  parsed <- parse_lsmeans_specs(specs)
  pairwise <- isTRUE(pairwise) || parsed$pairwise
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    cli::cli_abort("{.arg level} must be a number strictly between 0 and 1.")
  }
  if (!adjust %in% stats::p.adjust.methods) {
    cli::cli_abort("Unknown multiplicity adjustment {.val {adjust}}.")
  }

  information <- lsmeans_contrast_matrix(
    object,
    specs = parsed$variables,
    at = at
  )
  estimates <- lsmeans_estimates(object, information, level)
  if (!pairwise) {
    return(estimates)
  }

  contrasts <- lsmeans_pairs(
    object,
    information,
    specs = parsed$variables,
    level = level,
    adjust = adjust
  )
  structure(
    list(lsmeans = estimates, contrasts = contrasts),
    class = "lmm_lsmeans_list"
  )
}

#' @export
print.lmm_lsmeans_list <- function(x, ...) {
  cli::cli_h2("Estimated marginal means")
  print(x$lsmeans)
  cli::cli_h2("Pairwise contrasts")
  print(x$contrasts)
  invisible(x)
}
