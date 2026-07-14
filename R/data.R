#' Multi-site repeated-measures experiment
#'
#' An unbalanced experiment with three drugs in three centers and three
#' repeated measurement times. This is the multilocation repeated-measures
#' example described by Milliken and Johnson (2009, Section 28.3).
#'
#' @format A data frame with 153 rows and 5 variables:
#' \describe{
#'   \item{Center}{Center identifier with levels R, S, and T.}
#'   \item{Drug}{Drug identifier with levels 1, 2, and 3.}
#'   \item{Subject}{Subject identifier within each center and drug.}
#'   \item{Time}{Repeated-measurement occasion with levels 1, 2, and 3.}
#'   \item{Y}{Continuous response, including 28 missing values.}
#' }
#' @source Milliken, G. A., and Johnson, D. E. (2009). *Analysis of Messy Data,
#'   Volume 1: Designed Experiments* (2nd ed.), Section 28.3. Chapman and
#'   Hall/CRC. \doi{10.1201/EBK1584883340}
"multicentre"

#' Split-plot example from the SAS MIXED procedure documentation
#'
#' A balanced split-plot experiment with four blocks, three whole-plot levels,
#' and two subplot levels. SAS uses this data set in Example 79.1 of the
#' SAS/STAT 14.3 documentation for the MIXED procedure.
#'
#' @format A data frame with 24 rows and 4 variables:
#' \describe{
#'   \item{Block}{Block identifier with four levels.}
#'   \item{A}{Whole-plot treatment factor with three levels.}
#'   \item{B}{Subplot treatment factor with two levels.}
#'   \item{Y}{Continuous response.}
#' }
#' @source SAS Institute Inc. (2017). *SAS/STAT 14.3 User's Guide*,
#'   Example 79.1, Split-Plot Design.
#'   \url{https://support.sas.com/documentation/onlinedoc/stat/}
"sas_split_plot"

#' Repeated-measures growth example from the SAS MIXED documentation
#'
#' Growth measurements for 11 girls and 16 boys at ages 8, 10, 12, and 14.
#' SAS uses these Potthoff and Roy data in Example 79.2 of the SAS/STAT 14.3
#' documentation for the MIXED procedure.
#'
#' @format A data frame with 108 rows and 4 variables:
#' \describe{
#'   \item{Person}{Person identifier with 27 levels.}
#'   \item{Gender}{Gender factor with levels F and M.}
#'   \item{Age}{Age in years.}
#'   \item{y}{Continuous growth response.}
#' }
#' @source SAS Institute Inc. (2017). *SAS/STAT 14.3 User's Guide*,
#'   Example 79.2, Repeated Measures.
#'   \url{https://support.sas.com/documentation/onlinedoc/stat/}
"sas_growth"

#' Random-coefficients example from the SAS MIXED documentation
#'
#' Pharmaceutical stability measurements for three batches at six shelf ages.
#' SAS uses these simulated data in Example 79.5 of the SAS/STAT 14.3
#' documentation for the MIXED procedure.
#'
#' @format A data frame with 108 rows and 4 variables:
#' \describe{
#'   \item{Batch}{Batch identifier with three levels.}
#'   \item{Month}{Shelf age in months.}
#'   \item{Replicate}{Replicate identifier within a batch and month.}
#'   \item{Y}{Assay response, including 24 missing values.}
#' }
#' @source SAS Institute Inc. (2017). *SAS/STAT 14.3 User's Guide*,
#'   Example 79.5, Random Coefficients.
#'   \url{https://support.sas.com/documentation/onlinedoc/stat/}
"sas_random_coefficients"

#' Line-source irrigation example from the SAS MIXED documentation
#'
#' Winter-wheat yields for three cultivars in three blocks. Each plot contains
#' twelve subplots arranged on two sides of a line-source sprinkler. SAS uses
#' these data in Example 79.6 of the SAS/STAT 14.3 documentation to combine
#' three random-intercept terms with a four-band Toeplitz residual covariance.
#'
#' @format A data frame with 108 rows and 6 variables:
#' \describe{
#'   \item{Block}{Block identifier with three levels.}
#'   \item{Cult}{Winter-wheat cultivar with three levels.}
#'   \item{Sbplt}{Ordered subplot position from 1 to 12.}
#'   \item{Irrig}{Irrigation level with six levels.}
#'   \item{Dir}{Side of the sprinkler with levels North and South.}
#'   \item{Y}{Continuous yield response.}
#' }
#' @source SAS Institute Inc. (2017). *SAS/STAT 14.3 User's Guide*,
#'   Example 79.6, Line-Source Sprinkler Irrigation.
#'   \url{https://support.sas.com/documentation/onlinedoc/stat/}
"sas_line_source"
