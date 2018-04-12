##########################
# gumbel

# Description:
# Specifies the link function of gumbel distribution function in the contex of general, linear models
# Can be used in function glm

# Output
# Should be used in context with glm in the base package

gumbel <- function()
{
  linkfun <- function(mu) - log(- log(mu) )
  linkinv <- function(eta) exp(- exp(-eta))
  mu.eta <- function(eta) exp(- exp(-eta)) * exp(-eta)
  valideta <- function(eta) TRUE
  link <- paste("gumbel")
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta,
                 name = link), class = "link-glm")
}

#########################
# estSurv

# Description
# Estimates the discrete survival function given discrete hazard rates
# The hazard rates may or not depend on covariates, but the covariates have to be equal for all estimates of the hazard rate!
# Therefore the hazard rates only depend on time. 

# Input
# haz: Estimated hazard rates, which only depend on a time interval (numeric vector). 
# It is assumed that the estimates contain the intervals [0, a_1), [a_1, a_2), ..., [a_{q-1}, a_{q}) and the last interval is computed

# Output
# Survival function over all intervals [0, a_1), [a_1, a_2), ..., [a_{q-1}, a_{q}), [a_q, a_{\infty}

estSurv <- function (haz) {
  
  # Input checks
  if(!is.vector(haz)) {stop(
    "Argument *haz* is not a vector! Please specify a numeric vector of estimated hazard rates.")}
  if(!all(-sqrt(.Machine$double.eps) <= haz & 
          haz - 1 <= sqrt(.Machine$double.eps))) {
    stop("Argument *haz* is not a vector of probabilities! Please specify a numeric vector of estimated hazard rates.")}
  
  erg <- cumprod(1-haz)
  names(erg) <- paste("S(T=", 1:length(haz), ")", sep="")
  return(erg)
}

#########################
# estSurvCens

# Description
# Estimates the censoring survival function using life table estimator

estSurvCens <- function(dataSet, timeColumn, eventColumns) {

  # # Expand short data to long competing risks format
  # dataLongFormat  <- dataLongCompRisks (dataSet=dataSet, timeColumn=timeColumn, eventColumns=eventColumns, 
  #                    timeAsFactor=timeAsFactor)
  # 
  # # Construct censoring data
  # dataSetLongCens <- dataCensoring (dataSetLong=dataLongFormat, respColumn="y", 
  #                                   timeColumn="timeInt")

  # Convert to censoring format
  dataSetCens <- dataCensoringShort (dataSet=dataSet, eventColumns=eventColumns,
                                     timeColumn=timeColumn)

  # Estimate nonparametric survival function of censoring variable
  tempLifeTab <- lifeTable (dataSet=dataSetCens, timeColumn="timeCens", 
                            censColumn="yCens")
  preG <- tempLifeTab [[1]] [, "S"]
  GT <- c(1, preG)
  names(GT) <- paste("t=", 0:length(preG), sep="")
  return(GT)
}

#########################
# estMargProb

# Description
# Estimates the discrete marginal probability P(T=t) given discrete hazard rates
# The hazard rates may or not depend on covariates, but the covariates have to be equal for all estimates of the hazard rate!
# Therefore the hazard rates only depend on time. 

# Input
# haz: Estimated hazard rates, which only depend on a time interval (numeric vector). 
# It is assumed that the estimates contain the intervals [0, a_1), [a_1, a_2), ..., [a_{q-1}, a_{q}) and the last interval is computed

# Output
# Survival function over all intervals [0, a_1), [a_1, a_2), ..., [a_{q-1}, a_{q}), [a_q, a_{\infty}

estMargProb <- function (haz) {
  
  # Input checks
  if(!is.vector(haz)) {stop("Argument *haz* is not a vector! Please specify a numeric vector of estimated hazard rates.")}
  if(!all(-sqrt(.Machine$double.eps) <= haz & 
          haz - 1 <= sqrt(.Machine$double.eps))) {
    stop("Argument *haz* is not between 0 and 1! Please specify a numeric vector of probabilities.")
  }
  if(length(haz)==0) {stop("Argument *haz* is empty! Please specify a numeric vector of probabilities.")}
  
  if(length(haz)>1) {
    EstSurv <- estSurv(haz)
    EstProb <- c(haz[1], sapply(2:length(haz), function (x) haz [x] * EstSurv [x-1]))
    names(EstProb) <- paste("P(T=", 1:length(haz), ")", sep="")
  }
  else {
    EstProb <- 1
    names(EstProb) <- paste("P(T=", 1, ")", sep="")
  }

  return(EstProb)
}
