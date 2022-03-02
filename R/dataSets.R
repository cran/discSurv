#' Multiple Spell employment data
#' 
#' Subsample of 1000 persons from the national longitudinal survey of youth 1979 data. 
#' Included covariates are age, children, ethnicity, marital status and sex. 
#' The bivariate responses current state (spell) and discrete time interval (year) are 
#' the last two columns.
#' 
#' \itemize{
#'   \item Column "id" is defined as identification number for each person.
#' }
#' \itemize{
#'   \item Column "age" represents the time-varying age of each person in years.
#' }
#' \itemize{
#'   \item Column "child" consists of values \itemize{
#'   \item 0 - No children
#'   \item 1 - Individual has child/children
#'   }
#' }
#' \itemize{
#'   \item Column "ethnicity" consists of values \itemize{
#'   \item 1 - Hispanic
#'   \item 2 - Black
#'   \item 3 - Other
#'   }
#' }
#' \itemize{
#'   \item Column "marriage" consists of values \itemize{
#'   \item 1 - Never Married
#'   \item 2 - Currently married
#'   \item 3 - Other/Divorced
#'   }
#' }
#' \itemize{
#'   \item Column "sex" consists of values \itemize{
#'   \item 1 - Male
#'   \item 2 - Female
#'   }
#' }
#' \itemize{
#'   \item Column "spell" represents the time-varying employment status of each person. 
#'   Possible values are \itemize{
#'   \item 1 - Employed
#'   \item 2 - Unemployed
#'   \item 3 - Out of labor force
#'   \item 4 - In active forces
#'   \item 0 - Censored
#'   }
#' }
#' \itemize{
#'   \item Column "year" represents the discrete time intervals in years.
#' }
#' @name unempMultiSpell
#' @docType data
#' @author David Koehler \email{koehler@imbie.uni-bonn.de}
#' @source \href{https://www.nlsinfo.org/content/cohorts/NLSY97}{National Longitudinal Survey of Youth}
#' @usage data(unempMultiSpell)
#' @keywords data
NULL

#'Crash 2 competing risk data
#'
#'Adapted version of the crash2 trial data as availlable in the package Hmisc. Both death or survival and main cause of death are included. 
#'Death times are discretized into days. Included covariates are sex and age of patient, elapsed time between injury and 
#'hospitalization, type of injury, systolic blood pressure, heart rate, respiratory rate, central capillary 
#'refill time and total glascow coma score.
#'#' \itemize{
#'   \item Column "time" is time until death or hospital discharge/transfer in weeks.
#' }
#' \itemize{
#' \item Column "Status" denotes type of death\itemize{
#'   \item bleeding
#'   \item head injury
#'   \item vascular occlusion 
#'   \item multi organ failure
#'   \item other
#'   \item NA if patient is not dead (see "statusSE")
#'   }
#' }
#' \itemize{
#'   \item Column "statusSE" denotes death or discharge/transfer from hospital\itemize{
#'   \item 0 - Transfer/Discharge
#'   \item 1 - Death
#'   }
#' }
#' \itemize{
#'   \item Column "sex" denotes sex of the patient.
#' }
#' \itemize{
#'   \item Column "age" denotes age of the patient in years.
#' }
#' \itemize{
#'   \item Column "injurytime" gives time in hours between injury and hospitalization.
#' }
#' \itemize{
#'   \item Column "injurytype" denotes type of injury, one in \itemize{
#'   \item blunt
#'   \item penetrating
#'   \item blunt and penetrating
#'   }
#' }
#' \itemize{
#'   \item Column "sbp" denotes systolic blood pressure in mmHg.
#' }
#' \itemize{
#'   \item Column "rr" denotes respiratory rate per minute.
#' }
#' \itemize{
#'   \item Column "cc" denotes central capillary refill time in seconds.
#' }
#' \itemize{
#'   \item Column "hr" denotes heart rate per minute.
#' }
#' \itemize{
#'   \item Column "gcs" denotes total Glascow Coma Score.
#' }
#'@name crash2
#'@docType data
#'@author David Koehler \email{koehler@imbie.uni-bonn.de}
#'@source \link[Hmisc]{getHdata}
#'@usage data(crash2)
#'@keywords data
NULL
