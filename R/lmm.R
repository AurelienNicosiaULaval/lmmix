#' Control numerical optimization for `lmm()`
#'
#' @param optimizer Optimizer to use, either `"nlminb"` or `"optim"`.
#' @param optim_method Method passed to [stats::optim()] when `optimizer =
#'   "optim"`.
#' @param max_iter Maximum number of optimizer iterations.
#' @param rel_tol Relative convergence tolerance.
#' @param x_tol Parameter convergence tolerance for [stats::nlminb()].
#' @param initial Optional numeric vector of unconstrained starting values.
#' @param lower,upper Bounds for unconstrained covariance parameters.
#' @param deriv_method Numerical differentiation method passed to `numDeriv`.
#'
#' @return An object of class `lmm_control`.
#' @export
lmm_control <- function(
  optimizer = c("nlminb", "optim"),
  optim_method = "BFGS",
  max_iter = 1000L,
  rel_tol = 1e-8,
  x_tol = 1e-8,
  initial = NULL,
  lower = -20,
  upper = 20,
  deriv_method = "Richardson"
) {
  optimizer <- match.arg(optimizer)
  if (!is.numeric(max_iter) || length(max_iter) != 1L || max_iter < 1) {
    cli::cli_abort("{.arg max_iter} must be a positive number.")
  }
  if (!is.numeric(rel_tol) || length(rel_tol) != 1L || rel_tol <= 0) {
    cli::cli_abort("{.arg rel_tol} must be a positive number.")
  }
  if (!is.numeric(x_tol) || length(x_tol) != 1L || x_tol <= 0) {
    cli::cli_abort("{.arg x_tol} must be a positive number.")
  }
  if (!is.numeric(lower) || !is.numeric(upper) || any(lower >= upper)) {
    cli::cli_abort("{.arg lower} must be smaller than {.arg upper}.")
  }
  if (!is.null(initial) && !is.numeric(initial)) {
    cli::cli_abort("{.arg initial} must be NULL or a numeric vector.")
  }

  structure(
    list(
      optimizer = optimizer,
      optim_method = optim_method,
      max_iter = as.integer(max_iter),
      rel_tol = rel_tol,
      x_tol = x_tol,
      initial = initial,
      lower = lower,
      upper = upper,
      deriv_method = deriv_method
    ),
    class = "lmm_control"
  )
}

#' Fit a Gaussian linear mixed model
#'
#' `lmm()` explicitly evaluates a profiled ML or REML criterion for
#' `V = Z G Z' + R`. It can combine random effects with an independent,
#' compound-symmetric, AR(1), Toeplitz, or unstructured residual covariance.
#' Rows that are incomplete for any variable used by the fixed, random, or
#' repeated formula are removed before fitting.
#'
#' The `random` argument accepts one grouping formula. The `repeated` argument
#' accepts one ordering variable and one grouping expression. Each repeated
#' group may have at most one observation per ordering value. The current
#' likelihood assembles a dense marginal covariance matrix.
#'
#' @param data A data frame or tibble. Placing `data` first makes the function
#'   pipe-friendly.
#' @param formula A two-sided fixed-effects formula.
#' @param random Optional one-sided random-effects formula such as
#'   `~ 1 | subject` or `~ 1 + time | subject`.
#' @param repeated Optional one-sided repeated-measures formula such as
#'   `~ time | subject`.
#' @param structure Residual covariance structure: `"id"`, `"cs"`, `"ar1"`,
#'   `"toep"`, or `"un"`.
#' @param method Estimation method: restricted maximum likelihood (`"REML"`)
#'   or maximum likelihood (`"ML"`).
#' @param ddf Denominator degrees-of-freedom method. `"satterthwaite"` and
#'   `"residual"` are implemented. `"kenward-roger"` is reserved for future
#'   implementation.
#' @param control An object returned by `lmm_control()`.
#'
#' @return An object of class `lmm`. The object stores fixed and random-effect
#'   estimates, covariance components, fitted values, residuals, optimization
#'   diagnostics, the analysis data, and the model matrices.
#'
#' @references
#' Harville, D. A. (1977). Maximum likelihood approaches to variance component
#' estimation and to related problems. *Journal of the American Statistical
#' Association*, 72(358), 320-338.
#' \doi{10.1080/01621459.1977.10480998}
#'
#' LaMotte, L. R. (2007). A direct derivation of the REML likelihood function.
#' *Statistical Papers*, 48(2), 321-327.
#' \doi{10.1007/s00362-006-0335-6}
#'
#' Pinheiro, J. C., and Bates, D. M. (2000). *Mixed-Effects Models in S and
#' S-PLUS*. Springer. \doi{10.1007/b98882}
#'
#' @examples
#' fit <- lmm(
#'   data = nlme::Orthodont,
#'   formula = distance ~ age + Sex,
#'   random = ~ 1 | Subject
#' )
#' summary(fit)
#' @export
lmm <- function(
  data,
  formula,
  random = NULL,
  repeated = NULL,
  structure = c("id", "cs", "ar1", "toep", "un"),
  method = c("REML", "ML"),
  ddf = c("satterthwaite", "residual", "kenward-roger"),
  control = lmm_control()
) {
  call <- match.call()
  formula <- check_formula(formula, "formula")
  if (length(formula) != 3L) {
    cli::cli_abort("{.arg formula} must be a two-sided formula.")
  }
  if (!is.null(random)) {
    check_formula(random, "random")
  }
  if (!is.null(repeated)) {
    check_formula(repeated, "repeated")
  }
  if (!inherits(control, "lmm_control")) {
    cli::cli_abort("{.arg control} must be created by {.fn lmm_control}.")
  }

  structure <- match_choice(
    structure[[1L]],
    c("id", "cs", "ar1", "toep", "un"),
    "structure"
  )
  method <- toupper(match.arg(toupper(method[[1L]]), c("REML", "ML")))
  ddf <- match_choice(
    ddf[[1L]],
    c("satterthwaite", "residual", "kenward-roger"),
    "ddf"
  )
  if (ddf == "kenward-roger") {
    message <- paste(
      "Kenward-Roger inference is not implemented.",
      "Use {.val satterthwaite} or {.val residual}."
    )
    cli::cli_abort(message)
  }

  prepared <- prepare_analysis_data(data, formula, random, repeated)
  analysis_data <- prepared$data
  fixed_frame <- stats::model.frame(
    formula,
    data = analysis_data,
    na.action = stats::na.fail,
    drop.unused.levels = TRUE
  )
  fixed_terms <- stats::terms(fixed_frame)
  y <- stats::model.response(fixed_frame)
  if (!is.numeric(y) || is.matrix(y)) {
    cli::cli_abort("The response in {.arg formula} must be a numeric vector.")
  }
  x <- stats::model.matrix(fixed_terms, fixed_frame)
  if (qr(x)$rank < ncol(x)) {
    message <- paste(
      "The fixed-effects model matrix is rank deficient;",
      "simplify {.arg formula}."
    )
    cli::cli_abort(message)
  }
  if (nrow(x) <= ncol(x) && method == "REML") {
    cli::cli_abort("REML requires more observations than fixed-effect columns.")
  }

  random_design <- make_random_design(random, analysis_data)
  repeated_design <- make_repeated_design(repeated, analysis_data, structure)
  covariance_spec <- make_covariance_spec(
    random_design,
    repeated_design,
    structure
  )
  design <- list(
    x = x,
    y = as.numeric(y),
    random = random_design,
    repeated = repeated_design,
    covariance_spec = covariance_spec
  )

  optimization <- fit_covariance_parameters(design, method, control)
  gls <- gls_at_eta(optimization$eta, design)
  if (is.null(gls)) {
    cli::cli_abort("The fitted covariance matrix is not positive definite.")
  }
  random_effects <- estimate_blup(gls, design)
  marginal_fitted <- drop(x %*% gls$beta)
  conditional_fitted <- marginal_fitted
  if (!is.null(random_effects)) {
    conditional_fitted <- conditional_fitted + as.numeric(
      random_design$z %*% as.vector(t(random_effects))
    )
  }

  omitted <- if (length(prepared$omitted) > 0L) {
    stats::setNames(
      structure(prepared$omitted, class = "omit"),
      row.names(data)[prepared$omitted]
    )
  } else {
    NULL
  }

  result <- list(
    call = call,
    formula = formula,
    random_formula = random,
    repeated_formula = repeated,
    terms = fixed_terms,
    contrasts = attr(x, "contrasts"),
    xlevels = factor_xlevels(fixed_frame),
    model_frame = fixed_frame,
    data = analysis_data,
    original_row_index = prepared$row_index,
    na.action = omitted,
    design = design,
    coefficients = gls$beta,
    beta_vcov = gls$beta_vcov,
    eta = optimization$eta,
    eta_vcov = optimization$eta_vcov,
    hessian = optimization$hessian,
    covariance = gls$covariance,
    covariance_components = estimate_covariance_components(
      optimization$eta,
      optimization$eta_vcov,
      design
    ),
    random_effects = random_effects,
    marginal_fitted = marginal_fitted,
    fitted = conditional_fitted,
    residuals = as.numeric(y) - conditional_fitted,
    method = method,
    ddf = ddf,
    structure = structure,
    log_likelihood = -optimization$objective,
    convergence = list(
      code = optimization$code,
      message = optimization$message,
      optimizer = optimization$optimizer,
      iterations = optimization$iterations,
      evaluations = optimization$evaluations,
      hessian_positive_definite = optimization$hessian_positive_definite,
      hessian_eigenvalues = optimization$hessian_eigenvalues
    ),
    control = control
  )
  class(result) <- c("lmm", "lmmix")
  result
}
