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
  invisible(x)
}

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
  tibble::as_tibble(out)
}

#' @export
VarCorr.lmm <- function(object, ...) {
  object$covariance_components
}

#' @export
vcov.lmm <- function(object, ...) {
  object$beta_vcov
}

#' @export
logLik.lmm <- function(object, ...) {
  structure(
    object$log_likelihood,
    class = "logLik",
    df = length(object$coefficients) + length(object$eta),
    nobs = nobs(object)
  )
}

#' @export
AIC.lmm <- function(object, ..., k = 2) {
  -2 * as.numeric(logLik(object)) + k * attr(logLik(object), "df")
}

#' @export
BIC.lmm <- function(object, ...) {
  -2 * as.numeric(logLik(object)) +
    log(nobs(object)) * attr(logLik(object), "df")
}

#' @export
deviance.lmm <- function(object, ...) {
  -2 * as.numeric(logLik(object))
}

#' @export
nobs.lmm <- function(object, ...) {
  nrow(object$design$x)
}

#' @export
fitted.lmm <- function(object, type = c("conditional", "marginal"), ...) {
  type <- match.arg(type)
  if (type == "conditional") object$fitted else object$marginal_fitted
}

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

#' @export
model.matrix.lmm <- function(object, data = NULL, ...) {
  if (is.null(data)) object$design$x else model_matrix_newdata(object, data)
}

#' @export
model.frame.lmm <- function(formula, ...) {
  formula$model_frame
}

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

#' @export
terms.lmm <- function(x, ...) {
  x$terms
}

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
  colnames(out) <- paste0(
    format(c((1 - level) / 2, (1 + level) / 2) * 100, trim = TRUE),
    "%"
  )
  rownames(out) <- parm
  out
}
