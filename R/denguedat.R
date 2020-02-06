#' denguedat: Dengue fever reporting data from Puerto Rico
#'
#' Surveillance data from CDC Division of Vector-Borne Diseases.
#' 1990-2010 case reporting data included.
#' The first column, \code{onset_week}, indicates the week of symptom onset.
#' The second column, \code{report_week}, indicates the week of case report.
#'
#' @docType data
#'
#' @usage data(denguedat)
#'
#' @format A data frame.
#'
#' @keywords dengue
#'
#' @examples
#' data(denguedat)
#' nowcast <- NobBS(denguedat, as.Date("1990-04-09"),units="1 week",onset_date="onset_week",
#' report_date="report_week")
#' nowcast$estimates

"denguedat"
