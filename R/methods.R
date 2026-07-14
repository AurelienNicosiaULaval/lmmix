#' Extract fixed effects
#'
#' @param object A fitted model.
#' @param ... Additional arguments.
#' @export
fixef <- function(object, ...) {
  UseMethod("fixef")
}

#' Extract random effects
#'
#' @inheritParams fixef
#' @export
ranef <- function(object, ...) {
  UseMethod("ranef")
}

#' Extract covariance components
#'
#' @inheritParams fixef
#' @export
VarCorr <- function(object, ...) {
  UseMethod("VarCorr")
}

#' Standard methods for `lmm` fits
#'
#' Methods are provided for extracting likelihood quantities, fitted values,
#' residuals, model matrices and formulas, as well as for prediction, type III
#' testing and fixed-effect confidence intervals.
#'
#' Conditional fitted values and predictions include empirical random effects.
#' Marginal values use fixed effects only. In prediction from new data, known
#' grouping levels use their fitted random effect. New levels are rejected by
#' default and receive a random contribution of zero only when explicitly
#' allowed.
#'
#' @param x,object,formula An `lmm` object. The `formula` name follows the
#'   argument name of the `model.frame()` generic.
#' @param ... Additional fitted `lmm` models for likelihood-ratio comparison,
#'   or arguments passed to the corresponding method.
#' @param k Penalty per parameter used by `AIC()`.
#' @param type For `fitted()`, either `"conditional"` or `"marginal"`. For
#'   `residuals()`, one of `"response"`, `"pearson"`, or `"marginal"`. For
#'   `anova()`, only type `3` is implemented.
#' @param refit Whether REML models that differ in fixed effects are refitted
#'   with ML before comparison.
#' @param which Diagnostic plot to draw: standardized residuals against fitted
#'   values, a normal quantile-quantile plot, or observed against fitted values.
#' @param newdata Optional data frame used for prediction.
#' @param re.form `NULL` includes empirical random effects. Any non-`NULL`
#'   value requests a fixed-effects-only prediction.
#' @param allow.new.levels Whether conditional predictions may include new
#'   grouping levels. Their random contribution is zero when allowed.
#' @param na.action Missing-value action retained for prediction-method
#'   compatibility.
#' @param data Optional data frame used to construct a fixed-effects model
#'   matrix.
#' @param fixed.only Whether `formula()` returns only the fixed formula. If
#'   `FALSE`, it returns fixed, random, and repeated formulas in a list.
#' @param parm Fixed-effect parameters requested from `confint()`, supplied by
#'   name or position.
#' @param level Confidence level for fixed-effect intervals.
#' @param adjusted Whether `vcov()` returns the Kenward-Roger-adjusted
#'   covariance when available.
#'
#' @return The return value follows the corresponding base R generic.
#'
#' @name lmm-methods
NULL

#' @rdname lmm-methods
#' @export
print.lmm <- function(x, ...) {
  cli::cli_text("Linear mixed model fit by {x$method}")
  cli::cli_text("Formula: {.code {deparse(x$formula)}}")
  if (!is.null(x$random_formula)) {
    random_formulas <- normalize_random_formulas(x$random_formula)
    for (index in seq_along(random_formulas)) {
      label <- names(x$design$random$terms)[[index]]
      cli::cli_text(
        "Random ({label}): {.code {deparse(random_formulas[[index]])}}"
      )
    }
  }
  if (!is.null(x$repeated_formula)) {
    cli::cli_text(
      "Repeated: {.code {deparse(x$repeated_formula)}} ({toupper(x$structure)})"
    )
  } else {
    cli::cli_text("Residual structure: {toupper(x$structure)}")
  }
  cli::cli_text("Log-likelihood: {format(x$log_likelihood, digits = 6)}")
  cli::cli_text("Convergence code: {x$convergence$code}")
  invisible(x)
}

#' Summarize an `lmm` fit
#'
#' @param object An `lmm` object.
#' @param level Confidence level for fixed effects.
#' @param ... Additional arguments.
#'
#' @return An object of class `summary.lmm`.
#' @export
summary.lmm <- function(object, level = 0.95, ...) {
  result <- list(
    call = object$call,
    formula = object$formula,
    random_formula = object$random_formula,
    repeated_formula = object$repeated_formula,
    method = object$method,
    ddf = object$ddf,
    structure = object$structure,
    fixed = fixed_effects_table(object, level = level),
    covariance = object$covariance_components,
    type3 = type3_table(object),
    criteria = c(
      logLik = as.numeric(logLik(object)),
      AIC = AIC(object),
      BIC = BIC(object),
      deviance = deviance(object)
    ),
    convergence = object$convergence
  )
  class(result) <- "summary.lmm"
  result
}

#' @export
print.summary.lmm <- function(x, ...) {
  cli::cli_h1("Linear mixed model")
  cli::cli_text("Estimation: {x$method}")
  cli::cli_text("Denominator df: {x$ddf}")
  cli::cli_text("Residual covariance: {toupper(x$structure)}")
  cli::cli_h2("Fixed effects")
  print(x$fixed)
  cli::cli_h2("Type III tests")
  print(x$type3)
  cli::cli_h2("Covariance parameters")
  print(x$covariance)
  cli::cli_h2("Information criteria")
  print(x$criteria)
  if (x$convergence$code != 0L) {
    cli::cli_alert_warning("Optimizer convergence code: {x$convergence$code}")
  }
  if (!isTRUE(x$convergence$hessian_positive_definite)) {
    cli::cli_alert_warning("The likelihood Hessian is not positive definite.")
  }
  invisible(x)
}

#' Print an `lmmix` result table
#'
#' Printed column headings use human-readable labels without dots. The
#' underlying object retains its programmatic column names for compatibility
#' with `broom` and downstream R code.
#'
#' @param x An `lmmix` result table.
#' @param ... Additional arguments passed to the tibble print method.
#'
#' @return `x`, invisibly.
#' @keywords internal
#' @export
print.lmm_table <- function(x, ...) {
  out <- x
  names(out) <- lmm_display_names(names(out))
  table_classes <- intersect(class(out), c("tbl_df", "tbl", "data.frame"))
  class(out) <- if (length(table_classes) > 0L) {
    table_classes
  } else {
    "data.frame"
  }
  print(out, ...)
  invisible(x)
}

#' @rdname lmm-methods
#' @export
coef.lmm <- function(object, ...) {
  fixef(object)
}

#' @export
fixef.lmm <- function(object, ...) {
  object$coefficients
}

#' @export
ranef.lmm <- function(object, ...) {
  if (is.null(object$random_effects)) {
    return(as_lmm_table(tibble::tibble()))
  }

  out <- lapply(names(object$random_effects), function(label) {
    effects <- object$random_effects[[label]]
    table <- data.frame(
      group = rownames(effects),
      effects,
      row.names = NULL,
      check.names = FALSE
    )
    names(table)[[1L]] <- object$design$random$terms[[label]]$group_label
    as_lmm_table(tibble::as_tibble(table))
  })
  names(out) <- names(object$random_effects)
  if (length(out) == 1L) out[[1L]] else out
}

#' @export
VarCorr.lmm <- function(object, ...) {
  object$covariance_components
}

#' @rdname lmm-methods
#' @export
vcov.lmm <- function(object, adjusted = TRUE, ...) {
  if (isTRUE(adjusted)) object$beta_vcov else object$beta_vcov_model
}

#' @rdname lmm-methods
#' @export
logLik.lmm <- function(object, ...) {
  structure(
    object$log_likelihood,
    class = "logLik",
    df = length(object$coefficients) + length(object$eta),
    nobs = nobs(object)
  )
}

#' @rdname lmm-methods
#' @export
AIC.lmm <- function(object, ..., k = 2) {
  -2 * as.numeric(logLik(object)) + k * attr(logLik(object), "df")
}

#' @rdname lmm-methods
#' @export
BIC.lmm <- function(object, ...) {
  -2 * as.numeric(logLik(object)) +
    log(nobs(object)) * attr(logLik(object), "df")
}

#' @rdname lmm-methods
#' @export
deviance.lmm <- function(object, ...) {
  -2 * as.numeric(logLik(object))
}

#' @rdname lmm-methods
#' @export
nobs.lmm <- function(object, ...) {
  nrow(object$design$x)
}

#' @rdname lmm-methods
#' @export
fitted.lmm <- function(object, type = c("conditional", "marginal"), ...) {
  type <- match.arg(type)
  values <- if (type == "conditional") {
    object$fitted
  } else {
    object$marginal_fitted
  }
  restore_excluded(object, values)
}

#' @rdname lmm-methods
#' @export
residuals.lmm <- function(
  object,
  type = c("response", "pearson", "marginal"),
  ...
) {
  type <- match.arg(type)
  if (type == "marginal") {
    return(restore_excluded(
      object,
      object$design$y - object$marginal_fitted
    ))
  }
  if (type == "pearson") {
    return(restore_excluded(
      object,
      object$residuals / sqrt(diag(object$covariance$r))
    ))
  }
  restore_excluded(object, object$residuals)
}

#' @rdname lmm-methods
#' @export
plot.lmm <- function(
  x,
  which = c("residuals", "qq", "fitted"),
  ...
) {
  which <- match.arg(which)
  diagnostics <- tibble::tibble(
    .fitted = x$fitted,
    .std.resid = x$residuals / sqrt(diag(x$covariance$r)),
    .observed = x$design$y
  )

  if (which == "residuals") {
    return(
      ggplot2::ggplot(
        diagnostics,
        ggplot2::aes(x = .fitted, y = .std.resid)
      ) +
        ggplot2::geom_hline(
          yintercept = 0,
          linewidth = 0.4,
          linetype = 2,
          colour = "grey50"
        ) +
        ggplot2::geom_point(...) +
        ggplot2::labs(
          x = "Fitted values",
          y = "Standardized residuals",
          title = "Residuals versus fitted values"
        )
    )
  }

  if (which == "qq") {
    return(
      ggplot2::ggplot(
        diagnostics,
        ggplot2::aes(sample = .std.resid)
      ) +
        ggplot2::stat_qq(...) +
        ggplot2::stat_qq_line(colour = "grey50") +
        ggplot2::labs(
          x = "Theoretical quantiles",
          y = "Standardized residuals",
          title = "Normal Q-Q plot"
        )
    )
  }

  ggplot2::ggplot(
    diagnostics,
    ggplot2::aes(x = .fitted, y = .observed)
  ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linewidth = 0.4,
      linetype = 2,
      colour = "grey50"
    ) +
    ggplot2::geom_point(...) +
    ggplot2::labs(
      x = "Fitted values",
      y = "Observed values",
      title = "Observed versus fitted values"
    )
}

#' @rdname lmm-methods
#' @export
sigma.lmm <- function(object, ...) {
  sqrt(mean(diag(object$covariance$r)))
}

random_prediction <- function(object, newdata, allow.new.levels) {
  contributions <- lapply(names(object$design$random$terms), function(label) {
    random <- object$design$random$terms[[label]]
    parsed <- random$parsed
    group <- make_group(newdata, parsed$group_vars, "newdata")
    random_frame <- stats::model.frame(
      parsed$terms_formula,
      data = newdata,
      na.action = stats::na.pass
    )
    z_base <- stats::model.matrix(parsed$terms_formula, random_frame)
    fitted_effects <- object$random_effects[[label]]
    matched <- match(as.character(group), rownames(fitted_effects))
    unknown <- is.na(matched) & !is.na(group)
    if (any(unknown) && !isTRUE(allow.new.levels)) {
      levels <- unique(as.character(group[unknown]))
      cli::cli_abort(c(
        "New grouping level{?s} in random term {.field {label}}.",
        "i" = "Unknown level{?s}: {.val {levels}}.",
        "i" = paste(
          "Use {.arg allow.new.levels = TRUE} for population-level",
          "predictions."
        )
      ))
    }
    effects <- matrix(0, nrow = nrow(newdata), ncol = ncol(fitted_effects))
    known <- !is.na(matched)
    effects[known, ] <- fitted_effects[matched[known], , drop = FALSE]
    rowSums(z_base * effects)
  })
  Reduce(`+`, contributions, init = numeric(nrow(newdata)))
}

#' @rdname lmm-methods
#' @export
predict.lmm <- function(
  object,
  newdata = NULL,
  re.form = NULL,
  allow.new.levels = FALSE,
  na.action = stats::na.pass,
  ...
) {
  if (is.null(newdata)) {
    if (is.null(re.form)) fitted(object) else fitted(object, type = "marginal")
  } else {
    x <- model_matrix_newdata(object, newdata)
    prediction <- drop(x %*% object$coefficients)
    include_random <- is.null(re.form) && !is.null(object$random_effects)
    if (include_random) {
      prediction <- prediction + random_prediction(
        object,
        newdata,
        allow.new.levels = allow.new.levels
      )
    }
    prediction
  }
}

#' @rdname lmm-methods
#' @export
model.matrix.lmm <- function(object, data = NULL, ...) {
  if (is.null(data)) object$design$x else model_matrix_newdata(object, data)
}

#' @rdname lmm-methods
#' @export
model.frame.lmm <- function(formula, ...) {
  formula$model_frame
}

#' @rdname lmm-methods
#' @export
formula.lmm <- function(x, fixed.only = TRUE, ...) {
  if (isTRUE(fixed.only)) {
    x$formula
  } else {
    list(
      fixed = x$formula,
      random = x$random_formula,
      repeated = x$repeated_formula
    )
  }
}

#' @rdname lmm-methods
#' @export
terms.lmm <- function(x, ...) {
  x$terms
}

#' @rdname lmm-methods
#' @export
anova.lmm <- function(object, ..., type = 3, refit = TRUE) {
  dots <- list(...)
  if (length(dots) == 0L) {
    if (!identical(as.integer(type), 3L)) {
      cli::cli_abort("Only type III fixed-effect tests are implemented.")
    }
    out <- type3_table(object)
    attr(out, "ddf") <- object$ddf
    class(out) <- c("anova.lmm", class(out))
    return(out)
  }
  compare_lmm_models(c(list(object), dots), refit = refit)
}

formula_signature <- function(formula) {
  paste(deparse(formula), collapse = "")
}

random_signatures <- function(object) {
  formulas <- normalize_random_formulas(object$random_formula)
  if (is.null(formulas)) {
    character()
  } else {
    vapply(
      formulas,
      formula_signature,
      character(1L)
    )
  }
}

same_covariance_specification <- function(first, second) {
  identical(first$structure, second$structure) &&
    identical(
      formula_signature(first$repeated_formula),
      formula_signature(second$repeated_formula)
    ) &&
    identical(random_signatures(first), random_signatures(second))
}

fixed_model_nested <- function(smaller, larger, tolerance = 1e-8) {
  combined <- cbind(larger$design$x, smaller$design$x)
  combined_rank <- qr(combined, tol = tolerance)$rank
  larger_rank <- qr(larger$design$x, tol = tolerance)$rank
  combined_rank == larger_rank
}

residual_structure_nested <- function(smaller, larger) {
  if (identical(smaller, larger)) {
    return(TRUE)
  }
  allowed <- list(
    id = c("cs", "ar1", "toep", "un"),
    cs = c("toep", "un"),
    ar1 = c("toep", "un"),
    toep = "un",
    un = character()
  )
  larger %in% allowed[[smaller]]
}

random_model_nested <- function(smaller, larger) {
  small_terms <- smaller$design$random$terms %||% list()
  large_terms <- larger$design$random$terms %||% list()
  if (length(small_terms) == 0L) {
    return(TRUE)
  }
  all(vapply(small_terms, function(small_term) {
    candidates <- vapply(large_terms, function(large_term) {
      identical(small_term$parsed$group_vars, large_term$parsed$group_vars) &&
        all(small_term$term_names %in% large_term$term_names)
    }, logical(1L))
    any(candidates)
  }, logical(1L)))
}

covariance_model_nested <- function(smaller, larger) {
  same_repeated <- identical(
    formula_signature(smaller$repeated_formula),
    formula_signature(larger$repeated_formula)
  )
  same_repeated &&
    residual_structure_nested(smaller$structure, larger$structure) &&
    random_model_nested(smaller, larger)
}

refit_lmm_ml <- function(object) {
  action <- switch(object$na_action,
    na.omit = stats::na.omit,
    na.exclude = stats::na.exclude,
    na.fail = stats::na.fail
  )
  lmm(
    data = object$original_data,
    formula = object$formula,
    random = object$random_formula,
    repeated = object$repeated_formula,
    structure = object$structure,
    method = "ML",
    ddf = "satterthwaite",
    control = object$control,
    na.action = action
  )
}

validate_comparison_data <- function(models) {
  if (!all(vapply(models, inherits, logical(1L), what = "lmm"))) {
    message <- paste(
      "Every model supplied to {.fn anova} must inherit from",
      "{.cls lmm}."
    )
    cli::cli_abort(message)
  }
  response <- vapply(models, function(x) all.vars(x$formula)[[1L]], "")
  if (length(unique(response)) != 1L) {
    cli::cli_abort("Compared models must use the same response variable.")
  }
  reference_data <- models[[1L]]$original_data
  same_data <- vapply(models[-1L], function(x) {
    identical(x$original_data, reference_data) &&
      identical(x$original_row_index, models[[1L]]$original_row_index)
  }, logical(1L))
  if (!all(same_data)) {
    cli::cli_abort("Compared models must use the same data and observations.")
  }
}

compare_lmm_models <- function(models, refit) {
  validate_comparison_data(models)
  fixed_signatures <- vapply(
    models,
    function(x) formula_signature(x$formula),
    character(1L)
  )
  fixed_differ <- length(unique(fixed_signatures)) > 1L
  covariance_differ <- any(vapply(models[-1L], function(x) {
    !same_covariance_specification(models[[1L]], x)
  }, logical(1L)))
  if (fixed_differ && covariance_differ) {
    message <- paste(
      "Compare fixed effects and covariance structures in separate model",
      "sequences."
    )
    cli::cli_abort(message)
  }

  methods <- vapply(models, `[[`, character(1L), "method")
  refitted <- FALSE
  if (fixed_differ && any(methods == "REML")) {
    if (!isTRUE(refit)) {
      cli::cli_abort(
        "Models with different fixed effects must be compared with ML."
      )
    }
    cli::cli_inform("Refitting models with ML before comparison.")
    models <- lapply(models, refit_lmm_ml)
    refitted <- TRUE
  } else if (covariance_differ && any(methods != "ML")) {
    message <- paste(
      "Covariance models must be fitted with ML before",
      "comparison."
    )
    cli::cli_abort(message)
  }

  parameter_count <- vapply(
    models,
    function(x) attr(logLik(x), "df"),
    numeric(1L)
  )
  order <- order(parameter_count)
  models <- models[order]
  parameter_count <- parameter_count[order]
  fixed_signatures <- fixed_signatures[order]
  log_likelihood <- vapply(
    models,
    function(x) as.numeric(logLik(x)),
    numeric(1L)
  )
  chisq <- df <- p_value <- rep(NA_real_, length(models))

  for (index in 2:length(models)) {
    smaller <- models[[index - 1L]]
    larger <- models[[index]]
    nested <- if (fixed_differ) {
      same_covariance_specification(smaller, larger) &&
        fixed_model_nested(smaller, larger)
    } else {
      identical(fixed_signatures[[index - 1L]], fixed_signatures[[index]]) &&
        covariance_model_nested(smaller, larger)
    }
    if (!nested || parameter_count[[index]] <= parameter_count[[index - 1L]]) {
      cli::cli_abort("The supplied models are not strictly nested.")
    }
    chisq[[index]] <- 2 * (
      log_likelihood[[index]] - log_likelihood[[index - 1L]]
    )
    df[[index]] <- parameter_count[[index]] - parameter_count[[index - 1L]]
    p_value[[index]] <- stats::pchisq(
      max(chisq[[index]], 0),
      df = df[[index]],
      lower.tail = FALSE
    )
  }
  if (covariance_differ) {
    message <- paste(
      "Chi-squared reference distributions can be approximate for",
      "covariance parameters on the boundary."
    )
    cli::cli_warn(message)
  }

  out <- as_lmm_table(tibble::tibble(
    model = paste0("Model ", seq_along(models)),
    npar = parameter_count,
    AIC = vapply(models, AIC, numeric(1L)),
    BIC = vapply(models, BIC, numeric(1L)),
    logLik = log_likelihood,
    deviance = -2 * log_likelihood,
    Chisq = chisq,
    Df = df,
    p.value = p_value
  ))
  attr(out, "refitted") <- refitted
  class(out) <- c("anova.lmm_list", class(out))
  out
}

#' @export
print.anova.lmm <- function(x, ...) {
  label <- attr(x, "ddf") %||% "model"
  cli::cli_text("Type III tests with {label} degrees of freedom")
  NextMethod("print")
  invisible(x)
}

#' @rdname lmm-methods
#' @export
confint.lmm <- function(
  object,
  parm = names(object$coefficients),
  level = 0.95,
  ...
) {
  if (is.numeric(parm)) {
    parm <- names(object$coefficients)[parm]
  }
  unknown <- setdiff(parm, names(object$coefficients))
  if (length(unknown) > 0L) {
    cli::cli_abort("Unknown fixed-effect parameter{?s}: {.field {unknown}}.")
  }
  table <- fixed_effects_table(object, level = level)
  table <- table[match(parm, table$term), ]
  out <- as.matrix(table[c("conf.low", "conf.high")])
  colnames(out) <- c("Lower", "Upper")
  rownames(out) <- parm
  out
}
