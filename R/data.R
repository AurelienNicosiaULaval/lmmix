#' Multi-site repeated-measures experiment
#'
#' An unbalanced experiment with three drugs in three centers and three
#' repeated measurement times. The data reproduce Table 5.20 of the source
#' thesis and include the missing responses shown in that table.
#'
#' @format A data frame with 153 rows and 5 variables:
#' \describe{
#'   \item{Center}{Center identifier with levels R, S, and T.}
#'   \item{Drug}{Drug identifier with levels 1, 2, and 3.}
#'   \item{Subject}{Subject identifier within each center and drug.}
#'   \item{Time}{Repeated-measurement occasion with levels 1, 2, and 3.}
#'   \item{Y}{Continuous response, including 28 missing values.}
#' }
#' @source Mahsa Mohseni Bonab, *Programmation R et SAS pour modèles linéaires
#'   mixtes*, Table 5.20, p. 144, locally supplied final manuscript.
"multicentre"
