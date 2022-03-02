##############################################################################################
# Calculates life table estimator, Survival function and standard deviation without covariates
##############################################################################################

# Description
# Constructs life table and estimates hazard rate, survival function and corresponding standard errors

# Input
# dataShort: Original data set of format data.frame in short format
# timeColumn: Scalar character giving the column name of the time variable
# eventColumn: Scalar character giving the column name of the censoring variable

# Output
# data frame with columns
  # Rows: Each row of the table is an half open interval with width=1 (character)
  # n: Number of individuals at risk in a given time interval (integer)
  # events: Observed number of events in a given time interval (integer)
  # dropouts: Observed number of dropouts in a given time interval (integer)
  # atRisk: Estimated number of individuals at risk, corrected by dropouts
  # hazard: Estimated risk of death (without covariates) in a given time interval
  # seHazard: Estimated standard deviation of estimated hazard
  # S: Estimated survival curve
  # seS: Estimated standard deviation of estimated survival function



#' Life Table Construction and Estimates
#' 
#' Constructs a life table and estimates discrete hazards, survival
#' functions, discrete cumulative hazards and their standard errors without
#' covariates.
#' 
#' 
#' @param dataShort Original data in short format ("class data.frame").
#' @param timeColumn Name of the column with discrete survival times ("character vector").
#' @param eventColumn Gives the column name of the event indicator (1=observed,
#' 0=censored) ("character vector").
#' @param intervalLimits Optional names of the intervals for each row, e. g.
#' [a_0, a_1), [a_1, a_2), ..., [a_{q-1}, a_q) ("character vector")
#' @return List containing an object of class "data.frame" with following
#' columns \itemize{ \item{n} Number of individuals at risk in a given time
#' interval (integer) \item{events} Observed number of events in a given time
#' interval (integer) \item{dropouts} Observed number of dropouts in a given
#' time interval (integer) \item{atRisk} Estimated number of individuals at
#' risk, corrected by dropouts (numeric) \item{hazard} Estimated risk of death
#' (without covariates) in a given time interval \item{seHazard} Estimated
#' standard deviation of estimated hazard \item{S} Estimated survival curve
#' \item{seS} Estimated standard deviation of estimated survival function
#' \item{cumHazard} Estimated cumulative hazard function \item{seCumHazard}
#' Estimated standard deviation of the estimated cumulative hazard function
#' \item{margProb} Estimated marginal probability of event in time interval }
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' 
#' Matthias Schmid \email{matthias.schmid@@imbie.uni-bonn.de}
#' @references
#' \insertRef{tutzModelDisc}{discSurv} \cr\cr
#' \insertRef{lawlessLifetime}{discSurv}
#' @keywords survival
#' @examples
#' 
#' # Example with unemployment data
#' library(Ecdat)
#' data(UnempDur)
#' 
#' # Extract subset of all persons smaller or equal the median of age
#' UnempDurSubset <- subset(UnempDur, age <= median(UnempDur$age))
#' LifeTabUnempDur <- lifeTable(dataShort = UnempDurSubset, timeColumn = "spell", 
#' eventColumn = "censor1")
#' LifeTabUnempDur
#' 
#' # Example with monoclonal gammapothy data
#' library(survival)
#' head(mgus)
#' 
#' # Extract subset of mgus
#' subMgus <- mgus [mgus$futime<=median(mgus$futime), ]
#' 
#' # Transform time in days to intervals [0, 1), [1, 2), [2, 3), ... , [12460, 12461)
#' mgusInt <- subMgus
#' mgusInt$futime <- mgusInt$futime + 1
#' LifeTabGamma <- lifeTable(dataShort = mgusInt, timeColumn= "futime", eventColumn = "death")
#' head(LifeTabGamma$Output, 25)
#' plot(x = 1:dim(LifeTabGamma$Output)[1], y = LifeTabGamma$Output$hazard, type = "l", 
#' xlab = "Time interval", ylab = "Hazard", las = 1, 
#' main = "Life table estimated marginal discrete hazards")
#' 
#' @export lifeTable
lifeTable <- function(dataShort, timeColumn, eventColumn, intervalLimits = NULL) {
  
  # Input Checks
  if(!is.data.frame(dataShort)) {stop("Data set is not of type data.frame! Please convert the data!")}
  if(!(length(timeColumn) == 1 & is.character(timeColumn))) {stop("The column name is not correctly specified! The argument must be a scalar as character.")}
  if(!(length(eventColumn) == 1 & is.character(eventColumn))) {stop("The column name is not correctly specified! The argument must be a scalar as character.")}
  if(!(is.null (intervalLimits) )) {
    if(!(is.character(intervalLimits))) {stop("The interval borders are not in character format! Please give the appropriate format, e. g. [0, a_1), [a_1, a_2), ..., [a_{q-1}, a_q) with a_q beeing the number of intervals")}
    MaxTime <- max(dataShort [, timeColumn])
    if(!(length(intervalLimits) == MaxTime)) {stop("The interval borders have not the same length as the number of unique intervals! Please give the appropriate format, e. g. [0, a_1), [a_1, a_2), ..., [a_{q-1}, a_q) with a_q beeing the number of intervals")}
  }
  
  # Data transformation to life table
  atRiskInitial <- dim(dataShort) [1]
  dataShort <- dataShort [order(dataShort [, timeColumn]), ]
  formulaInput <- as.formula(paste(eventColumn, timeColumn, sep = "~"))
  events <- aggregate(formulaInput, data = dataShort, FUN = function (x) sum(x)) [, 2]
  dropouts <- aggregate(formulaInput, data = dataShort, FUN = function (x) sum(1-x)) [, 2]
  atRiskInput <- c(atRiskInitial, atRiskInitial - cumsum(events + dropouts))
  atRiskInput <- atRiskInput [-length(atRiskInput)]
  times <- as.numeric(names(table(as.numeric(as.character(dataShort [, timeColumn])))))
  Index <- which(diff(times) > 1)
  while(any(diff(times) > 1)) {
    Index <- which(diff(times) > 1) [1]
    atRiskInput <- c(atRiskInput[1:Index], 
                     atRiskInput[Index] - (events[Index] + dropouts[Index]), 
                     atRiskInput [(Index + 1):length(atRiskInput)])
    events <- c(events [1:Index], 0, events [(Index + 1):length(events)])
    dropouts <- c(dropouts [1:Index], 0, dropouts [(Index + 1):length(dropouts)])
    times <- c(times [1:Index], times[Index] + 1, times [(Index + 1):length(times)])
  }
  
  # Correction of first observed category:
  # It is assumed that no event and no dropouts occur in unobserved intervals
  # Each line must correspond to one interval
  if(times[1] != 1) {
    atRiskInput <- c(rep(atRiskInput [1], times[1] - 1), atRiskInput)
    events <- c(rep(0, times[1] - 1), events)
    dropouts <- c(rep(0, times[1] - 1), dropouts)
    times <- c(1:(times[1] - 1), times)
  }
  
  # Estimation of hazards, survival function cumulative hazard and standard deviations
  atRisk <- atRiskInput - dropouts / 2
  haz <- events / atRisk
  S <- cumprod(1 - haz)
  sehaz <- sqrt((haz - haz^2) / atRisk)
  seS <- S * sqrt(cumsum(haz / (1 - haz) / (atRisk)))
  cumHazard <- cumsum(haz)
  seCumHazard <- sqrt(cumsum(events / atRisk^2))
  
  # Construct interval borders
  if(is.null(intervalLimits)) {
    RowNames <- paste("[", c(0, times [-length(times)]), ", ", times, ")", sep = "")
  } else {
    RowNames <- intervalLimits
  }
  
  # Estimate marginal event probability P(T=t)
  margProb <- haz * c(1, S[-length(S)])

  # RowNames <- RowNames [-length(RowNames)]
  Output <- data.frame(n = atRiskInput, events = events, dropouts = dropouts, atRisk = atRisk,
             hazard = haz, seHazard = sehaz,
             S = S, seS = seS, 
             cumHazard = cumHazard, seCumHazard = seCumHazard, margProb = margProb,
             row.names = RowNames)
  
  # Exclude last row because theoretically the hazard will be 1 because all subjects die. 
  # However the estimator can be zero if no event occurs. Therefore it is not reliable
  # Output <- Output [-dim(Output) [1], ]
  Output <- list(Output = Output)
  class(Output) <- "discSurvLifeTable"
  return(Output)
}

#' @rdname lifeTable
#' @param x Object of class "discSurvLifeTable"("class discSurvLifeTable")
#' @param \dots Additional arguments to the print function
# #' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @keywords survival
#' @method print discSurvLifeTable
#' @export
print.discSurvLifeTable <- function(x, ...) {
  x <- x [[1]]
  for(i in 1:dim(x) [2]) {
    x [, i] <- round(x [, i], 4)
  }
  
  print(x)
}
