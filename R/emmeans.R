#' Support for `emmeans`
#'
#' These methods allow `emmeans::emmeans()` to recover the fixed-effects
#' reference grid and use model-specific Satterthwaite degrees of freedom.
#'
#' @param object An `lmm` object.
#' @param frame Model frame.
#' @param ... Additional arguments passed by `emmeans`.
#'
#' @return An `emmeans` data-recovery object.
#' @exportS3Method emmeans::recover_data
recover_data.lmm <- function(object, frame = object$model_frame, ...) {
  emmeans::recover_data(
    object$call,
    stats::delete.response(object$terms),
    object$na.action,
    frame = frame,
    ...
  )
}

#' @param trms,xlev,grid Arguments supplied by `emmeans`.
#'
#' @return A basis list used by `emmeans`.
#' @rdname recover_data.lmm
#' @exportS3Method emmeans::emm_basis
emm_basis.lmm <- function(object, trms, xlev, grid, ...) {
  model_frame <- stats::model.frame(
    trms,
    grid,
    na.action = stats::na.pass,
    xlev = xlev
  )
  x <- stats::model.matrix(
    trms,
    model_frame,
    contrasts.arg = object$contrasts
  )
  x <- x[, names(object$coefficients), drop = FALSE]

  list(
    X = x,
    bhat = unname(object$coefficients),
    nbasis = matrix(NA_real_, 1L, 1L),
    V = object$beta_vcov,
    dffun = function(k, dfargs) {
      dfargs$contrast_df(dfargs$object, k)
    },
    dfargs = list(object = object, contrast_df = contrast_df),
    misc = list()
  )
}
