#' Estimated Survival Function for Competing Risks
#' 
#' Computes the survival function S(T>t|x) based on estimated hazards of a competing risks model.
#' The discrete hazards may or may not depend on covariates. The covariates have to
#' be equal across all estimated hazards. Therefore the given discrete hazards
#' should only vary over time.
#' 
#' The argument \emph{hazards} must be given for all intervals [a_0, a_1), [a_1,
#' a_2), ..., [a_{q-1}, a_q), [a_q, Inf).
#' 
#' @param hazards Estimated discrete hazards ("numeric matrix"). 
#' Discrete hazards of each time interval are stored in the rows and the number of columns equal to the number of events.
#' @return Estimated survival probabilities ("numeric vector")
#' @note It is assumed that all time points up to the last interval [a_q, Inf)
#' are available. If not already present, these can be added manually.
#' @author Moritz Berger \email{moritz.berger@@imbie.uni-bonn.de} \cr \url{https://www.imbie.uni-bonn.de/personen/dr-moritz-berger/}
#' @seealso \code{\link{estSurv}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv}
#' @keywords competing_risks discrete_survival
#' @examples
#' 
#' # Example unemployment data
#' library(Ecdat)
#' data(UnempDur)
#' 
#' # Select subsample
#' subUnempDur <- UnempDur [1:100, ]
#' 
#' # Convert to long format
#' UnempLong <- dataLongCompRisks(dataShort = subUnempDur, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor4"))
#' head(UnempLong)
#' 
#' # Estimate continuation ratio model with logit link
#' vglmFit <- VGAM::vglm(formula = cbind(e0, e1, e2) ~ timeInt + age + logwage, data = UnempLong,
#' family = VGAM::multinomial(refLevel = "e0"))
#' 
#' # Estimate discrete survival function given age, logwage of first person
#' hazards <- VGAM::predictvglm(vglmFit, newdata = subset(UnempLong, obj == 1), type = "response")[,-1]
#' SurvivalFuncCondX <- estSurvCompRisks(rbind(hazards, 0.5))
#' SurvivalFuncCondX
#' 
#' 
#' @export estSurvCompRisks
estSurvCompRisks <- function (hazards) {
  
  # Input checks
  if(!is.matrix(hazards)) {stop(
    "Argument *hazards* is not a matrix! Please specify a matrix of estimated hazards with the time intervals in the rows 
    and the events in the columns.")}
  if(!all(-sqrt(.Machine$double.eps) <= hazards & 
          hazards - 1 <= sqrt(.Machine$double.eps))) {
    stop("Argument *hazards* is not a matrix of probabilities! Please specify a matrix of estimated hazards with the 
          time intervals in the rows and the events in the columns.")}
  
  overall_haz <- apply(hazards, 1, sum)
  erg <- cumprod(1 - overall_haz)
  names(erg) <- paste("S(T=", 1:nrow(hazards), ")", sep = "")
  return(erg)
}

#######################################################
#' Estimated Marginal Probabilities for Competing Risks
#' 
#' Estimates the marginal probability P(T = t,R = r|x) based on estimated discrete hazards of a competing risks model.
#' The discrete hazards may or may not depend on covariates. The covariates have to
#' be equal across all estimated discrete hazards. Therefore the given discrete hazards
#' should only vary over time.
#' 
#' The argument \emph{hazards} must be given for all intervals [a_0, a_1), [a_1,
#' a_2), ..., [a_{q-1}, a_q), [a_q, Inf).
#' 
#' @param hazards Estimated discrete hazards ("numeric matrix"). 
#' Discrete hazards of each time interval are stored in the rows and the number of columns equal to the number of events.
#' @return Estimated marginal probabilities ("numeric matrix")
#' @note It is assumed that all time points up to the last interval [a_q, Inf)
#' are available. If not already present, these can be added manually. 
#' In competing risk settings the marginal probabilities of the last theoretical interval 
#' depend on the assumptions on the discrete hazards in the last theoretical interval. 
#' However the estimation of all previous discrete intervals is not affected by those assumptions.
#' @author Moritz Berger \email{moritz.berger@@imbie.uni-bonn.de} \cr \url{https://www.imbie.uni-bonn.de/personen/dr-moritz-berger/}
#' @seealso \code{\link{estMargProb}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv}
#' @keywords competing risks discrete survival
#' @examples
#' 
#' # Example unemployment data
#' library(Ecdat)
#' data(UnempDur)
#' 
#' # Select subsample
#' subUnempDur <- UnempDur [1:100, ]
#' 
#' # Convert to long format
#' UnempLong <- dataLongCompRisks(dataShort = subUnempDur, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor4"))
#' head(UnempLong)
#' 
#' # Estimate continuation ratio model with logit link
#' vglmFit <- VGAM::vglm(formula = cbind(e0, e1, e2) ~ timeInt + age + logwage, data = UnempLong,
#' family = VGAM::multinomial(refLevel = "e0"))
#' 
#' # Estimate discrete survival function given age, logwage of first person
#' hazards <- VGAM::predictvglm(vglmFit, newdata = subset(UnempLong, obj == 1), type = "response")[,-1]
#' 
#' # Estimate marginal probabilities given age, logwage of first person
#' # Example 1
#' # Assumption: Discrete hazards in last theoretical interval are equal for both event types
#' MarginalProbCondX <- estMargProbCompRisks(rbind(hazards, 0.5))
#' MarginalProbCondX
#' all.equal(sum(MarginalProbCondX), 1) # TRUE: Marginal probabilities must sum to 1!
#' 
#' # Example 2
#' # Assumption: Discrete hazards in last theoretical interval are event1=, event2=
#' MarginalProbCondX2 <- estMargProbCompRisks(rbind(hazards, c(0.75, 0.25)))
#' MarginalProbCondX2
#' all.equal(sum(MarginalProbCondX2), 1) # TRUE: Marginal probabilities must sum to 1!
#' 
#' # Compare marginal probabilities given X
#' all.equal(MarginalProbCondX[1:5, ], MarginalProbCondX2[1:5, ])
#' all.equal(MarginalProbCondX[6, ], MarginalProbCondX2[6, ])
#' 
#' @export estMargProbCompRisks
estMargProbCompRisks <- function (hazards) {
  
  # Input checks
  if(!is.matrix(hazards)) {stop(
    "Argument *hazards* is not a matrix! Please specify a matrix of estimated hazards with the time intervals in the rows 
    and the events in the columns.")}
  if(!all(-sqrt(.Machine$double.eps) <= hazards & 
          hazards - 1 <= sqrt(.Machine$double.eps))) {
    stop("Argument *hazards* is not a matrix of probabilities! Please specify a matrix of estimated hazards with the 
         time intervals in the rows and the events in the columns.")}
  
  EstSurv <- estSurvCompRisks(hazards)
  EstProb <- matrix(NA, nrow=nrow(hazards), ncol = ncol(hazards))
  for(j in 1:ncol(hazards)){
    EstProb[,j] <- hazards[,j]*c(1,EstSurv[1:(nrow(hazards)-1)])
  }
  rownames(EstProb) <- paste("P(T=", 1:nrow(hazards), ")", sep = "")
  colnames(EstProb) <- colnames(hazards)
  
  return(EstProb)
}
