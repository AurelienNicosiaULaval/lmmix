#' Simulated multi-site repeated-measures experiment
#'
#' A simulated, unbalanced experiment with three drugs in three centers and
#' three repeated measurement times. The design follows the description in
#' Annex B of the package development brief. These are example data and are
#' not the original data used to calculate the published SAS reference values.
#'
#' @format A data frame with 162 rows and 5 variables:
#' \describe{
#'   \item{Center}{Center identifier with levels R, S, and T.}
#'   \item{Drug}{Drug identifier with levels 1, 2, and 3.}
#'   \item{Subject}{Subject identifier.}
#'   \item{Time}{Repeated-measurement occasion with levels 1, 2, and 3.}
#'   \item{Y}{Simulated continuous response, including eight missing values.}
#' }
#' @source Simulated by the package authors for documentation and testing.
"multicentre"
