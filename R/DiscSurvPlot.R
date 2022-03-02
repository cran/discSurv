###################################
#' Plot Estimated Survival Function 
#'
#' Generates a plot of an estimated survival function S(T>t|x) based on estimated discrete hazards.
#' 
#' @param hazards Estimated discrete hazards ("numeric vector")
#' @param ... Further arguments passed to \code{\link{plot}}.
#' @author Moritz Berger \email{moritz.berger@@imbie.uni-bonn.de} \cr \url{https://www.imbie.uni-bonn.de/personen/dr-moritz-berger/}
#' @seealso \code{\link{estSurv}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv}
#' @keywords discrete survival
#' @examples
#' # Example unemployment data
#' library(Ecdat)
#' data(UnempDur)
#' 
#' # Select subsample
#' subUnempDur <- UnempDur [1:100, ]
#' 
#' # Convert to long format
#' UnempLong <- dataLong(dataShort = subUnempDur, timeColumn = "spell", eventColumn = "censor1")
#' head(UnempLong)
#' 
#' # Estimate binomial model with logit link
#' Fit <- glm(formula = y ~ timeInt + age + logwage, data = UnempLong, family = binomial())
#' 
#' # Estimate discrete survival function given age, logwage of first person
#' Tmax   <- max(subUnempDur$spell)
#' UnempEval <- dataLong(dataShort = UnempDur[1,], timeColumn = "spell", eventColumn = "censor1", 
#' aggTimeFormat = TRUE, lastTheoInt = Tmax)
#' hazard <- predict(Fit, newdata = UnempEval, type = "response")
#'
#' plotSurv(hazard)
#'
#'
#'
#'@export plotSurv
plotSurv <- function (hazards, ...){
  est_surv <- estSurv(c(0, hazards))
  plot(0:(length(est_surv) - 1), est_surv, type = "s", ylab = "Estimated survival probability", xlab = "Time", 
       las = 1, ...)
}

###############################################
#' Plot Estimated Cumulative Incidence Function
#'
#' Generates a plot of an estimated cumulative incidence function P(T <= t, event=k | x) based on estimated hazards
#' of a discrete competing risks model or a discrete subdistribution hazard model.
#' 
#' @param hazards Numeric matrix (where each column represents one event) or vector of estimated hazards("numeric matrix").
#' @param eventFocus Column that represent the primary event ("integer vector"). 
#' Only applicable in the case of competing risks. 
#' @param ... Further arguments passed to \code{\link{plot}}.
#' @author Moritz Berger \email{moritz.berger@@imbie.uni-bonn.de} \cr \url{https://www.imbie.uni-bonn.de/personen/dr-moritz-berger/}
#' @seealso \code{\link{estSurv}}, \code{\link{estCumInz}}, \code{\link{compRisksGEE}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv}
#' @keywords discrete survival
#' @examples
#' 
#' # Example with unemployment data
#' library(Ecdat)
#' data(UnempDur)
#' 
#' # Select subsample
#' SubUnempDur <- UnempDur [1:100, ]
#' 
#' ################################
#' # Competing risks model 
#' 
#' # Estimate GEE models for all events
#' estGEE <- compRisksGEE(datShort = SubUnempDur, dataTransform = "dataLongCompRisks", 
#' corstr = "independence", formulaVariable =~ timeInt + age + ui + logwage * ui, 
#' eventColumns = c("censor1", "censor2", "censor3", "censor4"), timeColumn = "spell")
#' 
#' # Estimate hazards of all events given the covariates of third person
#' SubUnempDurLong <- dataLongCompRisks(dataShort = SubUnempDur, 
#' eventColumns = c("censor1", "censor2", "censor3", "censor4"), timeColumn = "spell") 
#' preds <- predict(estGEE, subset(SubUnempDurLong, obj == 3))
#' 
#' plotCumInc(preds, eventFocus = 3)
#' 
#' 
#' ###############################
#' # Subdistribution hazards model
#' 
#' # Convert to long format
#' SubUnempDurLong <- dataLongSubDist(dataShort = SubUnempDur, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor2", "censor3", "censor4"), eventFocus = "censor1")
#' 
#' # Estimate continuation ratio model with logit link
#' glmFit <- glm(formula = y ~ timeInt + age + ui + logwage * ui, data = SubUnempDurLong, 
#' family = binomial(), weights = SubUnempDurLong$subDistWeights)
#' 
#' # Estimated subdistribution hazard given the covariates of the third person
#' preds <- predict(glmFit, type = "response", newdata = subset(SubUnempDurLong, obj == 3))
#' 
#' plotCumInc(preds)
#' 
#'@export plotCumInc
plotCumInc <- function (hazards, eventFocus = NULL, ...){
  if(is.null(dim(hazards))){
    est_surv   <- estSurv(c(0, hazards))
    est_CumInc <- 1 - est_surv
  } else{
    est_CumInc <- estCumInz(hazards, eventFocus)
  }
  plot(0:length(est_CumInc), c(0, est_CumInc), type = "s", ylab = "Estimated cumulative incidence", xlab = "Time", 
       las = 1, ...)
}






