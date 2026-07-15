#' lmmix: Linear mixed models with correlated residuals
#'
#' `lmmix` fits Gaussian linear mixed models by explicitly evaluating a
#' profiled ML or REML criterion. It supports models that combine random
#' effects with correlated residual structures and provides Satterthwaite and
#' Kenward-Roger inference.
#'
#' @section Model scope:
#' The residual covariance may be independent, compound symmetric, AR(1),
#' full or fixed-band Toeplitz, or unstructured. A model may have no random
#' intercept, correlated random slopes, or multiple independent crossed or
#' nested grouping formulas. The
#' repeated-measures formula defines one ordering variable and the blocks of
#' the residual covariance.
#'
#' @section Inference and post-processing:
#' Fitted objects provide coefficient tests, type III tests, estimated
#' marginal means, pairwise contrasts, covariance-parameter intervals,
#' parametric-bootstrap likelihood-ratio tests, empirical BLUPs, predictions,
#' simulation, model updates, standard model methods, `ggplot2` diagnostic
#' plots, `broom` methods, and `emmeans` interoperability.
#'
#' @section Current limits:
#' The package fits univariate Gaussian responses. The likelihood factors
#' independent connected covariance components separately, but fitted objects
#' and some inference calculations retain dense marginal matrices. It is not a
#' general large-scale sparse mixed-model engine.
#'
#' @references
#' Patterson, H. D., and Thompson, R. (1971). Recovery of inter-block
#' information when block sizes are unequal. *Biometrika*, 58(3), 545-554.
#' \doi{10.1093/biomet/58.3.545}
#'
#' Harville, D. A. (1977). Maximum likelihood approaches to variance component
#' estimation and to related problems. *Journal of the American Statistical
#' Association*, 72(358), 320-338.
#' \doi{10.1080/01621459.1977.10480998}
#'
#' Satterthwaite, F. E. (1946). An approximate distribution of estimates of
#' variance components. *Biometrics Bulletin*, 2(6), 110-114.
#' \doi{10.2307/3002019}
#'
#' Fai, A. H.-T., and Cornelius, P. L. (1996). Approximate F-tests of multiple
#' degree of freedom hypotheses in generalized least squares analyses of
#' unbalanced split-plot experiments. *Journal of Statistical Computation and
#' Simulation*, 54(4), 363-378.
#' \doi{10.1080/00949659608811740}
#'
#' Kenward, M. G., and Roger, J. H. (1997). Small sample inference for fixed
#' effects from restricted maximum likelihood. *Biometrics*, 53(3), 983-997.
#' \doi{10.2307/2533558}
#'
#' Searle, S. R., Speed, F. M., and Milliken, G. A. (1980). Population marginal
#' means in the linear model: An alternative to least squares means.
#' *The American Statistician*, 34(4), 216-221.
#' \doi{10.1080/00031305.1980.10483031}
#'
#' @importFrom stats AIC BIC anova coef confint deviance fitted formula logLik
#' @importFrom stats model.frame model.matrix nobs predict residuals sigma
#' @importFrom stats simulate terms update vcov
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  ".fitted", ".observed", ".resid", ".std.resid", "component", "estimate",
  "group", "std.error", "term"
))
