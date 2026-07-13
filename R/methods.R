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
#' grouping levels use their fitted random effect and new levels receive a
#' random contribution of zero.
#'
#' @param x,object,formula An `lmm` object. The `formula` name follows the
#'   argument name of the `model.frame()` generic.
#' @param ... Additional arguments. Supplying another model to `anova()` is not
#'   supported.
#' @param k Penalty per parameter used by `AIC()`.
#' @param type For `fitted()`, either `"conditional"` or `"marginal"`. For
#'   `residuals()`, one of `"response"`, `"pearson"`, or `"marginal"`. For
#'   `anova()`, only type `3` is implemented.
#' @param which Diagnostic plot to draw: standardized residuals against fitted
#'   values, a normal quantile-quantile plot, or observed against fitted values.
#' @param newdata Optional data frame used for prediction.
#' @param re.form `NULL` includes empirical random effects. Any non-`NULL`
#'   value requests a fixed-effects-only prediction.
#' @param na.action Missing-value action retained for prediction-method
#'   compatibility.
#' @param data Optional data frame used to construct a fixed-effects model
#'   matrix.
#' @param fixed.only Whether `formula()` returns only the fixed formula. If
#'   `FALSE`, it returns fixed, random, and repeated formulas in a list.
#' @param parm Fixed-effect parameters requested from `confint()`, supplied by
#'   name or position.
#' @param level Confidence level for fixed-effect intervals.
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
    cli::cli_text("Random: {.code {deparse(x$random_formula)}}")
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
    return(tibble::tibble())
  }

  out <- data.frame(
    group = rownames(object$random_effects),
    object$random_effects,
    row.names = NULL,
    check.names = FALSE
  )
  names(out)[[1L]] <- object$design$random$group_label
  dot_names(tibble::as_tibble(out))
}

#' @export
VarCorr.lmm <- function(object, ...) {
  object$covariance_components
}

#' @rdname lmm-methods
#' @export
vcov.lmm <- function(object, ...) {
  object$beta_vcov
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
  if (type == "conditional") object$fitted else object$marginal_fitted
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
    return(object$design$y - object$marginal_fitted)
  }
  if (type == "pearson") {
    return(object$residuals / sqrt(diag(object$covariance$r)))
  }
  object$residuals
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
    .fitted = fitted(x),
    .std.resid = residuals(x, type = "pearson"),
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

random_prediction <- function(object, newdata) {
  random <- object$design$random
  parsed <- random$parsed
  group <- make_group(newdata, parsed$group_vars, "newdata")
  random_frame <- stats::model.frame(
    parsed$terms_formula,
    data = newdata,
    na.action = stats::na.pass
  )
  z_base <- stats::model.matrix(parsed$terms_formula, random_frame)
  matched <- match(as.character(group), rownames(object$random_effects))
  effects <- matrix(
    0,
    nrow = nrow(newdata),
    ncol = ncol(object$random_effects)
  )
  known <- !is.na(matched)
  effects[known, ] <- object$random_effects[matched[known], , drop = FALSE]
  rowSums(z_base * effects)
}

#' @rdname lmm-methods
#' @export
predict.lmm <- function(
  object,
  newdata = NULL,
  re.form = NULL,
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
      prediction <- prediction + random_prediction(object, newdata)
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
anova.lmm <- function(object, ..., type = 3) {
  dots <- list(...)
  if (length(dots) > 0L) {
    cli::cli_abort("Model-comparison ANOVA is not implemented.")
  }
  if (!identical(as.integer(type), 3L)) {
    cli::cli_abort("Only type III fixed-effect tests are implemented.")
  }
  out <- type3_table(object)
  attr(out, "ddf") <- object$ddf
  class(out) <- c("anova.lmm", class(out))
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
  colnames(out) <- c("conf.low", "conf.high")
  rownames(out) <- parm
  out
}
