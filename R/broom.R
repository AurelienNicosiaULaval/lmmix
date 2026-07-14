#' Tidy an `lmm` model
#'
#' @param x An `lmm` object.
#' @param effects Effects to return: fixed effects, covariance parameters, or
#'   both.
#' @param conf.int Whether to include confidence limits for fixed effects.
#' @param conf.level Confidence level.
#' @param ... Additional arguments.
#'
#' @return A tibble.
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.lmm <- function(
  x,
  effects = c("fixed", "ran_pars", "all"),
  conf.int = FALSE,
  conf.level = 0.95,
  ...
) {
  effects <- match.arg(effects)
  fixed <- fixed_effects_table(x, level = conf.level)
  if (!conf.int) {
    fixed$conf.low <- NULL
    fixed$conf.high <- NULL
  }
  fixed$effect <- "fixed"
  fixed <- fixed[c("effect", setdiff(names(fixed), "effect"))]
  if (effects == "fixed") {
    return(as_lmm_table(fixed))
  }

  covariance <- x$covariance_components
  covariance$term <- paste(covariance$group, covariance$term, sep = ":")
  covariance$effect <- "ran_pars"
  covariance$statistic <- NA_real_
  covariance$df <- NA_real_
  covariance$p.value <- NA_real_
  covariance$group <- NULL
  covariance$component <- NULL
  if (conf.int) {
    covariance$conf.low <- NA_real_
    covariance$conf.high <- NA_real_
  }
  covariance <- covariance[names(fixed)]

  if (effects == "ran_pars") {
    return(as_lmm_table(covariance))
  }
  as_lmm_table(tibble::as_tibble(rbind(fixed, covariance)))
}

#' One-row model summary for an `lmm` model
#'
#' @param x An `lmm` object.
#' @param ... Additional arguments.
#'
#' @return A one-row tibble.
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.lmm <- function(x, ...) {
  as_lmm_table(tibble::tibble(
    logLik = as.numeric(logLik(x)),
    AIC = AIC(x),
    BIC = BIC(x),
    deviance = deviance(x),
    df = attr(logLik(x), "df"),
    nobs = nobs(x),
    convergence = x$convergence$code,
    method = x$method
  ))
}

#' Augment data with fitted values and residuals
#'
#' @param x An `lmm` object.
#' @param data Data used for augmentation. Defaults to the analysis data.
#' @param newdata Optional new data. If supplied, it takes precedence over
#'   `data`.
#' @param ... Additional arguments passed to [predict()].
#'
#' @return A tibble with `.fitted`, `.resid`, and `.std.resid` columns.
#' @importFrom generics augment
#' @exportS3Method generics::augment
augment.lmm <- function(x, data = x$data, newdata = NULL, ...) {
  target <- newdata %||% data
  if (!is.data.frame(target)) {
    cli::cli_abort("Augmentation data must be a data frame or tibble.")
  }

  default_data <- is.null(newdata) && identical(target, x$data)
  fitted_values <- if (default_data) x$fitted else predict(x, target, ...)
  response_name <- all.vars(x$formula)[[1L]]
  residual_values <- if (response_name %in% names(target)) {
    target[[response_name]] - fitted_values
  } else {
    rep(NA_real_, nrow(target))
  }
  scale <- if (default_data) {
    sqrt(diag(x$covariance$r))
  } else {
    rep(sigma(x), nrow(target))
  }

  out <- tibble::as_tibble(target)
  out$.fitted <- fitted_values
  out$.resid <- residual_values
  out$.std.resid <- residual_values / scale
  as_lmm_table(out)
}
