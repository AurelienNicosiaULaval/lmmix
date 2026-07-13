#' lmmix: Linear mixed models with correlated residuals
#'
#' `lmmix` fits Gaussian linear mixed models by explicitly evaluating a
#' profiled ML or REML criterion. It supports models that combine random
#' effects with correlated residual structures and provides Satterthwaite
#' inference.
#'
#' @importFrom stats AIC BIC anova coef confint deviance fitted formula logLik
#' @importFrom stats model.frame model.matrix nobs predict residuals sigma terms
#' @importFrom stats vcov
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  ".fitted", ".resid", ".std.resid", "component", "estimate", "group",
  "std.error", "term"
))
