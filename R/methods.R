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
#' `ranef()` always returns a named list with one tibble per random-effect
#' term, including for models with a single term. It returns an empty list for
#' marginal models. `sigma()` is defined as `sqrt(mean(diag(R)))`; it is a
#' descriptive residual scale when residual variances differ across
#' observations.
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
#' @param test Reference test for likelihood-ratio model comparisons. The
#'   default `"chisq"` uses the usual asymptotic chi-squared reference. The
#'   `"parametric.bootstrap"` option simulates under each smaller model and is
#'   appropriate when covariance parameters can lie on a boundary.
#' @param nsim Number of simulations for a parametric-bootstrap comparison.
#' @param seed Optional integer seed for a reproducible parametric bootstrap.
#' @param which Diagnostic plot to draw: standardized residuals against fitted
#'   values, a normal quantile-quantile plot, or observed against fitted values.
#' @param newdata Optional data frame used for prediction.
#' @param re.form `NULL` includes every empirical random-effect term, `NA` or
#'   `~ 0` excludes all random effects, and a random-effects formula or list of
#'   formulas selects a subset of fitted terms.
#' @param allow.new.levels Whether conditional predictions may include new
#'   grouping levels. Their random contribution is zero when allowed.
#' @param se.fit Whether `predict()` returns standard errors for the expected
#'   response. These standard errors use the fixed-effect covariance and
#'   condition on any empirical random effects included through `re.form`.
#' @param interval Prediction interval type. `"confidence"` describes the
#'   expected response. `"prediction"` additionally includes the fitted
#'   residual variance but not uncertainty in empirical random effects.
#' @param na.action Missing-value action retained for prediction-method
#'   compatibility.
#' @param data Optional data frame used to construct a fixed-effects model
#'   matrix.
#' @param fixed.only Whether `formula()` returns only the fixed formula. If
#'   `FALSE`, it returns fixed, random, and repeated formulas in a list.
#' @param parm Parameters requested from `confint()`, supplied by name or fixed
#'   effect position. Use `"beta_"` for every fixed effect and `"theta_"` for
#'   every covariance parameter.
#' @param level Confidence level for fixed-effect and covariance intervals.
#' @param adjusted Whether `vcov()` returns the Kenward-Roger-adjusted
#'   covariance when available.
#' @param formula. Optional formula update accepted by `update()`.
#' @param evaluate Whether `update()` evaluates the updated call.
#'
#' @references
#' Self, S. G., and Liang, K.-Y. (1987). Asymptotic properties of maximum
#' likelihood estimators and likelihood ratio tests under nonstandard
#' conditions. *Journal of the American Statistical Association*, 82(398),
#' 605-610. \doi{10.1080/01621459.1987.10478472}
#'
#' Davison, A. C., and Hinkley, D. V. (1997). *Bootstrap Methods and Their
#' Application*. Cambridge University Press. \doi{10.1017/CBO9780511802843}
#'
#' @return The return value follows the corresponding base R generic.
#'
#' @name lmm-methods
NULL

#' @rdname lmm-methods
#' @export
print.lmm <- function(x, ...) {
  lines <- c(
    paste("Linear mixed model fit by", x$method),
    paste("Formula:", paste(deparse(x$formula), collapse = " "))
  )
  if (!is.null(x$random_formula)) {
    random_formulas <- normalize_random_formulas(x$random_formula)
    for (index in seq_along(random_formulas)) {
      label <- names(x$design$random$terms)[[index]]
      lines <- c(
        lines,
        paste0(
          "Random (", label, "): ",
          paste(deparse(random_formulas[[index]]), collapse = " ")
        )
      )
    }
  }
  if (!is.null(x$repeated_formula)) {
    lines <- c(
      lines,
      paste0(
        "Repeated: ", paste(deparse(x$repeated_formula), collapse = " "),
        " (", toupper(x$structure_label), ")"
      )
    )
  } else {
    lines <- c(
      lines,
      paste("Residual structure:", toupper(x$structure_label))
    )
  }
  lines <- c(
    lines,
    paste("Log-likelihood:", format(x$log_likelihood, digits = 6)),
    paste("Convergence code:", x$convergence$code)
  )
  writeLines(lines)
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
    structure = object$structure_label,
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
  writeLines(c(
    "Linear mixed model",
    paste("Estimation:", x$method),
    paste("Denominator df:", x$ddf),
    paste("Residual covariance:", toupper(x$structure)),
    "",
    "Fixed effects"
  ))
  print(x$fixed)
  writeLines("Type III tests")
  print(x$type3)
  writeLines("Covariance parameters")
  print(x$covariance)
  writeLines("Information criteria")
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
    return(list())
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
  out
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
      object$design$response - object$marginal_fitted
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
  rlang::check_installed(
    "ggplot2",
    reason = "to draw diagnostic plots for lmm models"
  )
  which <- match.arg(which)
  diagnostics <- tibble::tibble(
    .fitted = x$fitted,
    .std.resid = x$residuals / sqrt(diag(x$covariance$r)),
    .observed = x$design$response
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

selected_random_labels <- function(object, re.form) {
  available <- names(object$design$random$terms %||% list())
  if (is.null(re.form)) {
    return(available)
  }
  if (is.atomic(re.form) && length(re.form) == 1L && is.na(re.form)) {
    return(character())
  }
  if (inherits(re.form, "formula")) {
    right <- paste(deparse(re.form[[2L]]), collapse = "")
    if (right %in% c("0", "-1")) {
      return(character())
    }
  }
  if (length(available) == 0L) {
    cli::cli_abort("This model does not contain random-effect terms.")
  }

  requested <- normalize_random_formulas(re.form)
  fitted <- normalize_random_formulas(object$random_formula)
  requested_signatures <- vapply(requested, formula_signature, character(1L))
  fitted_signatures <- vapply(fitted, formula_signature, character(1L))
  matched <- match(requested_signatures, fitted_signatures)
  if (anyNA(matched)) {
    unknown <- requested_signatures[is.na(matched)]
    cli::cli_abort(
      paste0(
        "Unknown random-effect specification in {.arg re.form}: ",
        "{.code {unknown}}."
      )
    )
  }
  available[matched]
}

random_prediction <- function(object, newdata, allow.new.levels, labels) {
  contributions <- lapply(labels, function(label) {
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

fitted_random_prediction <- function(object, labels) {
  contributions <- lapply(labels, function(label) {
    random <- object$design$random$terms[[label]]
    effects <- as.vector(t(object$random_effects[[label]]))
    as.numeric(random$z %*% effects)
  })
  Reduce(`+`, contributions, init = numeric(nrow(object$design$x)))
}

prediction_residual_variance <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    return(diag(object$covariance$r))
  }
  base <- object$covariance$residual_base
  if (object$structure == "id") {
    variance <- rep(base[[1L]], nrow(newdata))
  } else {
    time_variable <- object$design$repeated$time_variable
    if (is.null(time_variable) || !time_variable %in% names(newdata)) {
      cli::cli_abort(
        "Prediction intervals require the repeated-measures ordering variable."
      )
    }
    time <- match(
      as.character(newdata[[time_variable]]),
      object$design$repeated$time_levels
    )
    if (anyNA(time)) {
      cli::cli_abort(
        "Prediction intervals cannot use unseen repeated-measures time values."
      )
    }
    variance <- diag(base)[time]
  }
  weights <- prediction_weights(object, newdata)
  variance / weights
}

prediction_offset <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    return(object$offset)
  }
  fixed_terms <- stats::delete.response(object$terms)
  frame <- stats::model.frame(
    fixed_terms,
    data = newdata,
    xlev = object$xlevels,
    na.action = stats::na.pass
  )
  formula_offset <- stats::model.offset(frame)
  if (is.null(formula_offset)) {
    formula_offset <- numeric(nrow(newdata))
  }
  explicit <- if (is.null(object$offset_expression)) {
    numeric(nrow(newdata))
  } else {
    evaluate_data_vector(
      object$offset_expression,
      newdata,
      object$evaluation_environment,
      "offset"
    )
  }
  as.numeric(formula_offset) + explicit
}

prediction_weights <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    return(object$weights)
  }
  if (is.null(object$weights_expression)) {
    return(rep(1, nrow(newdata)))
  }
  weights <- evaluate_data_vector(
    object$weights_expression,
    newdata,
    object$evaluation_environment,
    "weights"
  )
  if (any(!is.finite(weights) | weights <= 0)) {
    cli::cli_abort("{.arg weights} must contain positive finite values.")
  }
  weights
}

#' @rdname lmm-methods
#' @export
predict.lmm <- function(
  object,
  newdata = NULL,
  re.form = NULL,
  allow.new.levels = FALSE,
  se.fit = FALSE,
  interval = c("none", "confidence", "prediction"),
  level = 0.95,
  na.action = stats::na.pass,
  ...
) {
  interval <- match.arg(interval)
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit)) {
    cli::cli_abort("{.arg se.fit} must be TRUE or FALSE.")
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    cli::cli_abort("{.arg level} must be one number strictly between 0 and 1.")
  }
  labels <- selected_random_labels(object, re.form)
  if (is.null(newdata)) {
    x <- object$design$x
    prediction <- object$marginal_fitted +
      fitted_random_prediction(object, labels)
  } else {
    x <- model_matrix_newdata(object, newdata)
    prediction <- prediction_offset(object, newdata) +
      drop(x %*% object$coefficients)
    if (length(labels) > 0L) {
      prediction <- prediction + random_prediction(
        object,
        newdata,
        allow.new.levels = allow.new.levels,
        labels = labels
      )
    }
  }

  if (!isTRUE(se.fit) && interval == "none") {
    return(prediction)
  }
  fixed_variance <- rowSums((x %*% object$beta_vcov) * x)
  fixed_se <- sqrt(pmax(fixed_variance, 0))
  df <- vapply(seq_len(nrow(x)), function(index) {
    contrast_df(object, x[index, ])
  }, numeric(1L))
  fit <- prediction
  if (interval != "none") {
    interval_variance <- fixed_variance
    if (interval == "prediction") {
      interval_variance <- interval_variance + prediction_residual_variance(
        object,
        newdata = newdata
      )
    }
    critical <- stats::qt((1 + level) / 2, df = df)
    width <- critical * sqrt(pmax(interval_variance, 0))
    fit <- cbind(
      fit = prediction,
      lwr = prediction - width,
      upr = prediction + width
    )
  }
  if (isTRUE(se.fit)) {
    return(list(fit = fit, se.fit = fixed_se, df = df))
  }
  fit
}

#' @rdname lmm-methods
#' @export
simulate.lmm <- function(object, nsim = 1, seed = NULL, ...) {
  if (
    !is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
      nsim < 1 || nsim > .Machine$integer.max || nsim != floor(nsim)
  ) {
    cli::cli_abort("{.arg nsim} must be one positive integer.")
  }
  restore_seed <- set_bootstrap_seed(seed)
  if (!is.null(restore_seed)) {
    on.exit(restore_seed(), add = TRUE)
  }
  chol_v <- chol_factor(object$covariance$v)
  if (is.null(chol_v)) {
    cli::cli_abort("The fitted covariance matrix is not positive definite.")
  }
  nsim <- as.integer(nsim)
  simulations <- replicate(nsim, {
    object$marginal_fitted +
      drop(t(chol_v) %*% stats::rnorm(nobs(object)))
  })
  simulations <- as.data.frame(simulations, optional = TRUE)
  names(simulations) <- paste0("sim_", seq_len(nsim))
  simulations
}

#' @rdname lmm-methods
#' @export
update.lmm <- function(object, formula., ..., evaluate = TRUE) {
  formula <- if (missing(formula.)) {
    object$formula
  } else {
    stats::update.formula(object$formula, formula.)
  }
  arguments <- list(
    data = object$original_data,
    formula = formula,
    random = object$random_formula,
    repeated = object$repeated_formula,
    structure = object$structure_label,
    method = object$method,
    ddf = object$ddf,
    control = object$control,
    na.action = switch(object$na_action,
      na.omit = stats::na.omit,
      na.exclude = stats::na.exclude,
      na.fail = stats::na.fail
    )
  )
  if (!is.null(object$weights_input)) {
    arguments$weights <- object$weights_input
  }
  if (!is.null(object$offset_input)) {
    arguments$offset <- object$offset_input
  }
  if (!is.null(object$contrasts_input)) {
    arguments$contrasts <- object$contrasts_input
  }
  changes <- list(...)
  if (length(changes) > 0L) {
    arguments[names(changes)] <- changes
  }
  call <- as.call(c(list(quote(lmm)), arguments))
  if (!isTRUE(evaluate)) {
    return(call)
  }
  do.call(lmm, arguments)
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
anova.lmm <- function(
  object,
  ...,
  type = 3,
  refit = TRUE,
  test = c("chisq", "parametric.bootstrap"),
  nsim = 199L,
  seed = NULL
) {
  dots <- list(...)
  test <- match.arg(test)
  if (length(dots) == 0L) {
    if (test != "chisq") {
      cli::cli_abort(
        "{.arg test} applies only to comparisons of two or more models."
      )
    }
    if (!identical(as.integer(type), 3L)) {
      cli::cli_abort("Only type III fixed-effect tests are implemented.")
    }
    out <- type3_table(object)
    attr(out, "ddf") <- object$ddf
    class(out) <- c("anova.lmm", class(out))
    return(out)
  }
  compare_lmm_models(
    c(list(object), dots),
    refit = refit,
    test = test,
    nsim = nsim,
    seed = seed
  )
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
    identical(first$covariance_order, second$covariance_order) &&
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
  small_structure <- smaller$structure
  large_structure <- larger$structure
  small_spec <- smaller$design$covariance_spec
  large_spec <- larger$design$covariance_spec

  if (identical(small_structure, large_structure)) {
    if (small_structure != "toep") {
      return(TRUE)
    }
    return(small_spec$covariance_order <= large_spec$covariance_order)
  }
  if (small_structure == "id") {
    return(large_structure %in% c("cs", "ar1", "toep", "un"))
  }
  if (small_structure %in% c("cs", "ar1")) {
    if (large_structure == "un") {
      return(TRUE)
    }
    return(
      large_structure == "toep" &&
        large_spec$covariance_order == large_spec$n_times
    )
  }
  small_structure == "toep" && large_structure == "un"
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
    residual_structure_nested(smaller, larger) &&
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
    structure = object$structure_label,
    method = "ML",
    ddf = "satterthwaite",
    weights = object$weights_input,
    offset = object$offset_input,
    contrasts = object$contrasts_input,
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

bootstrap_lrt_fit <- function(object, response) {
  design <- object$design
  design$y <- response
  control <- object$control
  control$initial <- object$eta
  tryCatch(
    suppressWarnings(fit_covariance_parameters(
      design,
      method = "ML",
      control = control,
      compute_hessian = FALSE
    )),
    error = function(cnd) NULL
  )
}

bootstrap_lrt_pair <- function(smaller, larger, observed, nsim) {
  mean <- drop(smaller$design$x %*% smaller$coefficients)
  chol_v <- chol_factor(smaller$covariance$v)
  if (is.null(chol_v)) {
    cli::cli_abort("The null-model covariance matrix is not positive definite.")
  }

  statistics <- rep(NA_real_, nsim)
  for (simulation in seq_len(nsim)) {
    response <- mean + drop(t(chol_v) %*% stats::rnorm(length(mean)))
    small_fit <- bootstrap_lrt_fit(smaller, response)
    large_fit <- bootstrap_lrt_fit(larger, response)
    converged <- !is.null(small_fit) && !is.null(large_fit) &&
      small_fit$code == 0L && large_fit$code == 0L &&
      is.finite(small_fit$objective) && is.finite(large_fit$objective)
    if (converged) {
      statistics[[simulation]] <- max(
        2 * (small_fit$objective - large_fit$objective),
        0
      )
    }
  }

  successful <- sum(is.finite(statistics))
  if (successful == 0L) {
    cli::cli_abort("Every parametric-bootstrap refit failed.")
  }
  if (successful < nsim) {
    cli::cli_warn(
      "{nsim - successful} of {nsim} parametric-bootstrap refits failed."
    )
  }
  valid <- statistics[is.finite(statistics)]
  list(
    p.value = (1 + sum(valid >= observed)) / (successful + 1),
    statistics = valid,
    successful = successful
  )
}

set_bootstrap_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }
  if (
    !is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed)
  ) {
    cli::cli_abort("{.arg seed} must be NULL or one non-negative integer.")
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  set.seed(as.integer(seed))
  function() {
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }
}

compare_lmm_models <- function(models, refit, test, nsim, seed) {
  validate_comparison_data(models)
  if (test == "parametric.bootstrap") {
    if (
      !is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
        nsim < 1 || nsim > .Machine$integer.max || nsim != floor(nsim)
    ) {
      cli::cli_abort("{.arg nsim} must be one positive integer.")
    }
    restore_seed <- set_bootstrap_seed(seed)
    if (!is.null(restore_seed)) {
      on.exit(restore_seed(), add = TRUE)
    }
    nsim <- as.integer(nsim)
  }
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
  bootstrap <- vector("list", length(models))
  if (test == "parametric.bootstrap") {
    for (index in 2:length(models)) {
      bootstrap[[index]] <- bootstrap_lrt_pair(
        models[[index - 1L]],
        models[[index]],
        observed = max(chisq[[index]], 0),
        nsim = nsim
      )
      p_value[[index]] <- bootstrap[[index]]$p.value
    }
  } else if (covariance_differ) {
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
  attr(out, "test") <- test
  if (test == "parametric.bootstrap") {
    attr(out, "nsim") <- nsim
    attr(out, "bootstrap") <- bootstrap
  }
  class(out) <- c("anova.lmm_list", class(out))
  out
}

#' @export
print.anova.lmm_list <- function(x, ...) {
  if (identical(attr(x, "test"), "parametric.bootstrap")) {
    writeLines(
      paste(
        "Likelihood-ratio tests with",
        attr(x, "nsim"),
        "bootstrap simulations"
      )
    )
  }
  NextMethod("print")
  invisible(x)
}

#' @export
print.anova.lmm <- function(x, ...) {
  label <- attr(x, "ddf") %||% "model"
  writeLines(paste("Type III tests with", label, "degrees of freedom"))
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
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    cli::cli_abort("{.arg level} must be one number strictly between 0 and 1.")
  }
  fixed_names <- names(object$coefficients)
  covariance <- covariance_natural(object$eta, object$design)
  covariance_names <- names(covariance)
  if (is.numeric(parm)) {
    parm <- fixed_names[parm]
  }
  parm <- unlist(lapply(parm, function(value) {
    switch(value,
      beta_ = fixed_names,
      theta_ = covariance_names,
      value
    )
  }), use.names = FALSE)
  parm <- unique(parm)
  known <- c(fixed_names, covariance_names)
  unknown <- setdiff(parm, known)
  if (length(unknown) > 0L) {
    cli::cli_abort("Unknown parameter{?s}: {.field {unknown}}.")
  }

  fixed_table <- fixed_effects_table(object, level = level)
  fixed_intervals <- as.matrix(fixed_table[c("conf.low", "conf.high")])
  rownames(fixed_intervals) <- fixed_table$term

  covariance_intervals <- matrix(
    NA_real_,
    nrow = length(covariance),
    ncol = 2L,
    dimnames = list(covariance_names, c("conf.low", "conf.high"))
  )
  standard_error <- object$covariance_components$std.error
  critical <- stats::qnorm(1 - (1 - level) / 2)
  components <- object$covariance_components$component
  for (index in seq_along(covariance)) {
    estimate <- unname(covariance[[index]])
    se <- standard_error[[index]]
    if (!is.finite(estimate) || !is.finite(se)) {
      next
    }
    if (components[[index]] == "var" && estimate > 0) {
      se_link <- se / estimate
      covariance_intervals[index, ] <- exp(
        log(estimate) + c(-1, 1) * critical * se_link
      )
    } else if (components[[index]] == "cor" && abs(estimate) < 1) {
      se_link <- se / (1 - estimate^2)
      covariance_intervals[index, ] <- tanh(
        atanh(estimate) + c(-1, 1) * critical * se_link
      )
    }
  }

  intervals <- rbind(fixed_intervals, covariance_intervals)
  out <- intervals[match(parm, rownames(intervals)), , drop = FALSE]
  colnames(out) <- c("Lower", "Upper")
  rownames(out) <- parm
  out
}
