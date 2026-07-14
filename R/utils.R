`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

as_lmm_table <- function(x) {
  class(x) <- unique(c("lmm_table", class(x)))
  x
}

lmm_display_names <- function(x) {
  labels <- c(
    effect = "Effect",
    term = "Term",
    estimate = "Estimate",
    std.error = "Std Error",
    statistic = "Statistic",
    df = "DF",
    p.value = "p value",
    conf.low = "Conf Low",
    conf.high = "Conf High",
    num.df = "Num DF",
    den.df = "Den DF",
    group = "Group",
    component = "Component",
    contrast = "Contrast",
    .fitted = "Fitted",
    .resid = "Residual",
    .std.resid = "Std Residual",
    logLik = "Log Lik",
    deviance = "Deviance",
    nobs = "N Obs",
    convergence = "Convergence",
    method = "Method"
  )
  display <- trimws(gsub("\\.+", " ", x))
  matched <- x %in% names(labels)
  display[matched] <- unname(labels[x[matched]])
  display
}

check_formula <- function(x, arg) {
  if (!inherits(x, "formula")) {
    cli::cli_abort("{.arg {arg}} must be a formula.")
  }
  x
}

parse_bar_formula <- function(x, arg) {
  check_formula(x, arg)
  if (length(x) != 2L) {
    cli::cli_abort("{.arg {arg}} must be a one-sided formula.")
  }

  bar <- x[[2L]]
  if (!is.call(bar) || !identical(bar[[1L]], as.name("|"))) {
    cli::cli_abort(
      "{.arg {arg}} must use the form {.code ~ terms | group}."
    )
  }

  group_vars <- all.vars(bar[[3L]])
  if (length(group_vars) == 0L) {
    cli::cli_abort("{.arg {arg}} must contain a grouping variable.")
  }

  list(
    terms_expr = bar[[2L]],
    terms_formula = rlang::new_formula(
      lhs = NULL,
      rhs = bar[[2L]],
      env = environment(x)
    ),
    group_expr = bar[[3L]],
    group_vars = group_vars
  )
}

make_group <- function(data, vars, label) {
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0L) {
    cli::cli_abort(
      "Unknown variable{?s} in {.arg {label}}: {.field {missing_vars}}."
    )
  }

  interaction(
    data[vars],
    drop = TRUE,
    lex.order = TRUE,
    sep = ":"
  )
}

required_variables <- function(formula, random, repeated) {
  unique(c(
    all.vars(formula),
    if (!is.null(random)) all.vars(random),
    if (!is.null(repeated)) all.vars(repeated)
  ))
}

prepare_analysis_data <- function(data, formula, random, repeated) {
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame or tibble.")
  }
  if (nrow(data) == 0L) {
    cli::cli_abort("{.arg data} must contain at least one row.")
  }

  vars <- required_variables(formula, random, repeated)
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0L) {
    cli::cli_abort("Unknown variable{?s}: {.field {missing_vars}}.")
  }

  used <- stats::complete.cases(data[vars])
  if (!any(used)) {
    cli::cli_abort(
      "No complete observations remain after removing missing data."
    )
  }

  list(
    data = data[used, , drop = FALSE],
    row_index = which(used),
    omitted = which(!used),
    variables = vars
  )
}

factor_xlevels <- function(model_frame) {
  factors <- vapply(model_frame, is.factor, logical(1L))
  lapply(model_frame[factors], levels)
}

chol_factor <- function(x) {
  tryCatch(chol(x), error = function(cnd) NULL)
}

chol_solve <- function(chol_x, rhs) {
  rhs <- as.matrix(rhs)
  backsolve(chol_x, forwardsolve(t(chol_x), rhs))
}

chol_logdet <- function(chol_x) {
  2 * sum(log(diag(chol_x)))
}

symmetrize <- function(x) {
  (x + t(x)) / 2
}

generalized_inverse <- function(x, tolerance = sqrt(.Machine$double.eps)) {
  eig <- eigen(symmetrize(x), symmetric = TRUE)
  cutoff <- tolerance * max(1, max(abs(eig$values)))
  keep <- eig$values > cutoff
  if (!any(keep)) {
    return(matrix(NA_real_, nrow(x), ncol(x)))
  }

  inverse <- eig$vectors[, keep, drop = FALSE] %*%
    diag(1 / eig$values[keep], nrow = sum(keep)) %*%
    t(eig$vectors[, keep, drop = FALSE])
  symmetrize(inverse)
}

invert_hessian <- function(hessian) {
  eig <- eigen(symmetrize(hessian), symmetric = TRUE, only.values = TRUE)$values
  tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(eig)))
  positive_definite <- all(is.finite(eig)) && all(eig > tolerance)

  list(
    vcov = generalized_inverse(hessian),
    positive_definite = positive_definite,
    eigenvalues = eig
  )
}

match_choice <- function(x, choices, arg) {
  if (length(x) != 1L || is.na(x)) {
    cli::cli_abort("{.arg {arg}} must have length one.")
  }
  match.arg(tolower(x), choices)
}

model_matrix_newdata <- function(object, newdata) {
  fixed_terms <- stats::delete.response(object$terms)
  model_frame <- stats::model.frame(
    fixed_terms,
    data = newdata,
    xlev = object$xlevels,
    na.action = stats::na.pass
  )
  stats::model.matrix(
    fixed_terms,
    data = model_frame,
    contrasts.arg = object$contrasts
  )
}
