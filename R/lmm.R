#' Control numerical optimization for `lmm()`
#'
#' @param optimizer Optimizer strategy: `"auto"`, `"nlminb"`, or `"optim"`.
#' @param optim_method Method passed to [stats::optim()] when `optimizer =
#'   "optim"`.
#' @param max_iter Maximum number of optimizer iterations.
#' @param rel_tol Relative convergence tolerance.
#' @param x_tol Parameter convergence tolerance for [stats::nlminb()].
#' @param initial Optional numeric vector of unconstrained starting values.
#' @param lower,upper Bounds for unconstrained covariance parameters.
#' @param deriv_method Numerical differentiation method passed to `numDeriv`.
#' @param max_restarts Maximum number of deterministic restart attempts after
#'   an unsuccessful initial fit.
#' @param restart_scale Size of deterministic perturbations to starting values.
#'
#' @return An object of class `lmm_control`.
#' @export
lmm_control <- function(
  optimizer = c("auto", "nlminb", "optim"),
  optim_method = "BFGS",
  max_iter = 1000L,
  rel_tol = 1e-8,
  x_tol = 1e-8,
  initial = NULL,
  lower = -20,
  upper = 20,
  deriv_method = "Richardson",
  max_restarts = 3L,
  restart_scale = 0.25
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
  if (
    !is.numeric(max_restarts) ||
      length(max_restarts) != 1L ||
      max_restarts < 0
  ) {
    cli::cli_abort("{.arg max_restarts} must be a non-negative number.")
  }
  if (
    !is.numeric(restart_scale) ||
      length(restart_scale) != 1L ||
      restart_scale < 0
  ) {
    cli::cli_abort("{.arg restart_scale} must be a non-negative number.")
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
      deriv_method = deriv_method,
      max_restarts = as.integer(max_restarts),
      restart_scale = restart_scale
    ),
    class = "lmm_control"
  )
}

#' Fit a Gaussian linear mixed model
#'
#' `lmm()` explicitly evaluates a profiled ML or REML criterion for
#' `V = Z G Z' + R`. It can combine random effects with an independent,
#' compound-symmetric, AR(1), full or fixed-band Toeplitz, or unstructured
#' residual covariance.
#' Missing values are handled according to `na.action`.
#'
#' The `random` argument accepts a grouping formula or a list of independent
#' grouping formulas. The `repeated` argument
#' accepts one ordering variable and one grouping expression. Each repeated
#' group may have at most one observation per ordering value. The current
#' likelihood assembles a dense marginal covariance matrix.
#'
#' @param data A data frame or tibble. Placing `data` first makes the function
#'   pipe-friendly.
#' @param formula A two-sided fixed-effects formula.
#' @param random Optional one-sided random-effects formula such as
#'   `~ 1 | subject`, or a list of such formulas.
#' @param repeated Optional one-sided repeated-measures formula such as
#'   `~ time | subject`.
#' @param structure Residual covariance structure: `"id"`, `"cs"`, `"ar1"`,
#'   `"toep"`, `"toep(k)"`, or `"un"`. The `"toep"` form estimates every
#'   available Toeplitz band. The `"toep(k)"` form estimates `k` bands,
#'   including the main diagonal, and fixes longer-lag covariances to zero.
#' @param method Estimation method: restricted maximum likelihood (`"REML"`)
#'   or maximum likelihood (`"ML"`).
#' @param ddf Denominator degrees-of-freedom method: `"satterthwaite"`,
#'   `"kenward-roger"`, or `"residual"`.
#' @param control An object returned by `lmm_control()`.
#' @param na.action Missing-value action: [stats::na.omit()],
#'   [stats::na.exclude()], or [stats::na.fail()].
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
#' Kenward, M. G., and Roger, J. H. (1997). Small sample inference for fixed
#' effects from restricted maximum likelihood. *Biometrics*, 53(3), 983-997.
#' \doi{10.2307/2533558}
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
  control = lmm_control(),
  na.action = stats::na.omit
) {
  call <- match.call()
  formula <- check_formula(formula, "formula")
  if (length(formula) != 3L) {
    cli::cli_abort("{.arg formula} must be a two-sided formula.")
  }
  random_formulas <- normalize_random_formulas(random)
  if (!is.null(random_formulas)) {
    lapply(random_formulas, check_formula, arg = "random")
  }
  if (!is.null(repeated)) {
    check_formula(repeated, "repeated")
  }
  if (!inherits(control, "lmm_control")) {
    cli::cli_abort("{.arg control} must be created by {.fn lmm_control}.")
  }

  structure_spec <- parse_covariance_structure(structure[[1L]])
  structure <- structure_spec$name
  method <- toupper(match.arg(toupper(method[[1L]]), c("REML", "ML")))
  ddf <- match_choice(
    ddf[[1L]],
    c("satterthwaite", "residual", "kenward-roger"),
    "ddf"
  )
  if (ddf == "kenward-roger" && method != "REML") {
    cli::cli_abort("Kenward-Roger inference requires {.val REML} estimation.")
  }

  prepared <- prepare_analysis_data(
    data,
    formula,
    random_formulas,
    repeated,
    na.action = na.action
  )
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

  random_design <- make_random_design(random_formulas, analysis_data)
  repeated_design <- make_repeated_design(repeated, analysis_data, structure)
  covariance_spec <- make_covariance_spec(
    random_design,
    repeated_design,
    structure,
    covariance_order = structure_spec$order
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
  kr <- NULL
  beta_vcov <- gls$beta_vcov
  if (ddf == "kenward-roger") {
    kr <- kenward_roger_adjustment(
      optimization$eta,
      optimization$eta_vcov,
      design,
      gls$beta_vcov
    )
    beta_vcov <- kr$adjusted_vcov
  }
  marginal_fitted <- drop(x %*% gls$beta)
  conditional_fitted <- marginal_fitted
  if (!is.null(random_effects)) {
    for (label in names(random_design$terms)) {
      conditional_fitted <- conditional_fitted + as.numeric(
        random_design$terms[[label]]$z %*%
          as.vector(t(random_effects[[label]]))
      )
    }
  }

  omitted <- if (length(prepared$omitted) > 0L) {
    action_class <- if (prepared$na_action == "na.exclude") {
      "exclude"
    } else {
      "omit"
    }
    stats::setNames(
      structure(prepared$omitted, class = action_class),
      row.names(data)[prepared$omitted]
    )
  } else {
    NULL
  }

  result <- list(
    call = call,
    formula = formula,
    random_formula = if (inherits(random, "formula")) {
      random
    } else {
      random_formulas
    },
    repeated_formula = repeated,
    terms = fixed_terms,
    contrasts = attr(x, "contrasts"),
    xlevels = factor_xlevels(fixed_frame),
    model_frame = fixed_frame,
    data = analysis_data,
    original_data = data,
    original_row_index = prepared$row_index,
    na.action = omitted,
    na_action = prepared$na_action,
    design = design,
    coefficients = gls$beta,
    beta_vcov = beta_vcov,
    beta_vcov_model = gls$beta_vcov,
    kr = kr,
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
    covariance_order = covariance_spec$covariance_order,
    structure_label = covariance_spec$structure_label,
    log_likelihood = -optimization$objective,
    convergence = list(
      code = optimization$code,
      message = optimization$message,
      optimizer = optimization$optimizer,
      iterations = optimization$iterations,
      evaluations = optimization$evaluations,
      hessian_positive_definite = optimization$hessian_positive_definite,
      hessian_eigenvalues = optimization$hessian_eigenvalues,
      attempts = optimization$attempts
    ),
    control = control
  )
  class(result) <- c("lmm", "lmmix")
  result
}
