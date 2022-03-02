########################################
# GEE model for discrete competing risks

#' GEE model for discrete competing risks
#' 
#' Estimates generalized estimation equation model for each competing event separately. 
#' Dependence within person IDs is accounted for by assuming a working covariance structure.
#' 
#' @param datShort Original data set in short format with each row corresponding to one independent 
#' observation("class data.frame").
#' @param dataTransform Specification of the data transformation function from short to long format("character vector"). 
#' There are two available options: Without time dependent covariates ("dataLongCompRisks") and with 
#' time dependent covariates ("dataLongCompRisksTimeDep"). The default is set to the former.
#' @param corstr Assumption of correlation structure ("character vector"). The following are 
#' permitted: '"independence"', '"exchangeable"', '"ar1"', '"unstructured"' and '"userdefined".
#' @param formulaVariable Specifies the right hand side of the regression formula ("class formula").
#' The default is to use the discrete time variable, which corresponds to a covariate free hazard. 
#' It is recommended to always include the discrete time variable "timeInt".
#' @param \dots Additional arguments to data transformation (compRisksGEE) or prediction function (predict).
#' Preprocessing function argument responseAsFactor has to be set to FALSE (Default).
#' @return Returns an object of class "geeglm".
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @details Variables in argument \emph{formulaVariable} need to be separated by "+ ". 
#' For example if the two variables \emph{timeInt} and \emph{X1} should be included the formula would be
#' "~ timeInt + X1". The variable \emph{timeInt} is constructed before estimation of the model.
#' @seealso \code{\link{covarGEE}}, \code{\link{dataLongCompRisks}}, \code{\link{dataLongCompRisksTimeDep}}, 
#' \code{\link[geepack]{geeglm}}
#' @references 
#' \insertRef{minjungDiscComp}{discSurv}
#' @keywords survival
#' @examples
#' 
#' # Example with unemployment data
#' library(Ecdat)
#' data(UnempDur)
#' 
#' # Select subsample
#' SubUnempDur <- UnempDur [1:100, ]
#' 
#' # Estimate GEE models for all events
#' estGEE <- compRisksGEE(datShort = SubUnempDur, dataTransform = "dataLongCompRisks", 
#' corstr = "independence", formulaVariable =~ timeInt + age + ui + logwage * ui, 
#' eventColumns = c("censor1", "censor2", "censor3", "censor4"), timeColumn = "spell")
#' names(estGEE)
#' estGEE[[1]]
#' 
#' # Predictions
#' SubUnempDurLong <- dataLongCompRisks(dataShort = SubUnempDur, 
#' eventColumns = c("censor1", "censor2", "censor3", "censor4"), timeColumn = "spell")
#' preds <- predict(estGEE, newdata = SubUnempDurLong)
#' head(preds)
#' 
#' @export compRisksGEE
compRisksGEE <- function(datShort, dataTransform = "dataLongCompRisks", corstr = "independence", 
                             formulaVariable =~ timeInt, ... ){

  # Data transformation
  if(dataTransform=="dataLongCompRisks"){
    datLong <- dataLongCompRisks(dataShort = datShort, responseAsFactor=FALSE, ... )
  }
  if(dataTransform == "dataLongCompRisksTimeDep"){
    datLong <- dataLongCompRisksTimeDep(dataSemiLong = datShort, responseAsFactor = FALSE, ... )
  }

  # Construct formula
  eventColumns <- grep("e[0-9]+", names(datLong), value = TRUE)
  # Delete censoring column
  eventColumns <- eventColumns[-which(eventColumns == "e0")]
  eventColumns_length <- length(eventColumns)
  constructFormula <- sapply(1:eventColumns_length, 
          function(x) formula(paste(paste(eventColumns[x], collapse = ","), 
          "~", paste(attr(terms(formulaVariable),"term.labels"), collapse=" + "), sep = " ")))

  # Fit GEE
  RES <- vector(mode = "list", length = eventColumns_length)
  names(RES) <- paste(eventColumns)
  for(k in 1:eventColumns_length){
    RES[[k]] <-  geeglm(formula = constructFormula[[k]], family = binomial, data = datLong, 
                        corstr = corstr, id = datLong$obj)
  }
  class(RES) <- "dCRGEE"
  return(RES)
}

###############################################################################

#' @rdname compRisksGEE
#' @param object Discrete time competing risks GEE model prediction model ("class dCRGEE").
# #' @author Thomas Welchowski
#' @param newdata ("class data.set") New data set to be used for prediction (class data.frame).
#' @keywords survival
#' @method predict dCRGEE
#' @export
predict.dCRGEE <- function (object, newdata, ...) {

  predMat <- matrix(NA, nrow = dim(newdata)[1], ncol = length(object))
  for(k in 1:length(object)){
    
    predMat[, k] <- predict(object[[k]], newdata=newdata, type = "response", ...)
    
  }
  dimnames(predMat)[[2]] <- names(object)
  return(predMat)
}

################################################################################
# Internal GEE functions

#' @rdname compRisksGEE
#' @param k Competing cause k ("integer vector")
#' @param i Cluster number i ("integer vector"), e. g. patients or hospitals
#' @param l Discrete time interval l ("integer vector")
#' @param modelEst Discrete time competing risks GEE model prediction model ("class dCRGEE").
# #' @author Thomas Welchowski
#' @keywords survival
#' @noRd
gradientLogLikGEE <- function(k, i, l, modelEst){
  
  currentModel <- modelEst[[k]]
  
  exp_eta <- exp(currentModel$linear.predictors[, 1])
  modelMat <- model.matrix(currentModel$formula, currentModel$model)
  dimModelMat <- dim(modelMat)
  
  gradAllObs <- matrix(NA, nrow=dimModelMat[1], ncol=dimModelMat[2])
  for( j in 1:dim(modelMat)[2] ){
    
    varInterest <- modelMat[, j]
    
    gradAllObs[, j] <- currentModel$y * varInterest / ( 1 + exp_eta ) - 
      (1-currentModel$y) * varInterest * exp_eta / ( 1 + exp_eta )
    
    
  }
  
  dimnames(gradAllObs)[[2]] <- dimnames(modelMat)[[2]]
  gradAllObsSub <- gradAllObs[currentModel$data$obj==i & currentModel$data$timeInt==l, ]
  
  gradAllObsSub
}

#' @rdname compRisksGEE
#' @param k Competing event k ("integer vector")
#' @param i Cluster number i ("integer vector"), e. g. patients or hospitals
#' @param modelEst Discrete time competing risks GEE model prediction model ("class dCRGEE").
# #' @author Thomas Welchowski
#' @keywords survival
#' @noRd
U_k_i <- function(k, i, modelEst){
  
  currentModel <- modelEst[[k]]
  modelMat <- model.matrix(currentModel$formula, currentModel$model)
  n_i <- max(currentModel$data[currentModel$data$obj==i, "timeInt"])
  rowSums(sapply(1:n_i, function(m) gradientLogLikGEE(k=k, i=i, l=m, modelEst=modelEst)))
  
}

#' @rdname compRisksGEE
#' @param k Competing event k ("integer vector")
#' @param i Cluster number i ("integer vector"), e. g. patients or hospitals
#' @param l Discrete time interval l ("integer vector")
#' @param modelEst Discrete time competing risks GEE model prediction model ("class dCRGEE").
# #' @author Thomas Welchowski
#' @keywords survival
#' @noRd
hessianLogitLikGEE <- function(k, i, l, modelEst){
  
  currentModel <- modelEst[[k]]
  
  exp_eta <- exp(currentModel$linear.predictors[, 1])
  modelMat <- model.matrix(currentModel$formula, currentModel$model)
  dimModelMat <- dim(modelMat)
  
  restrictObs <- currentModel$data$obj==i & currentModel$data$timeInt==l
  hessMat <- matrix(NA, nrow=dimModelMat[2], ncol=dimModelMat[2])
  for( j in 1:dim(modelMat)[2] ){
    for( m in 1:dim(modelMat)[2] ){
      
      varInterest1 <- modelMat[restrictObs, j]
      varInterest2 <- modelMat[restrictObs, m]
      
      hessMat[m, j] <- - varInterest1 * varInterest2 * exp_eta[restrictObs] / 
        ( 1 + exp_eta[restrictObs] )^2
      
    }
  }
  
  dimnames(hessMat)[[1]] <- dimnames(modelMat)[[2]]
  dimnames(hessMat)[[2]] <- dimnames(modelMat)[[2]]
  hessMat
  
}

#' @rdname compRisksGEE
#' @param k Competing event k ("integer vector")
#' @param modelEst Discrete time competing risks GEE model prediction model ("class dCRGEE").
# #' @author Thomas Welchowski
#' @keywords survival
#' @noRd
H_k <- function(k, modelEst){
  
  currentModel <- modelEst[[k]]
  modelMat <- model.matrix(currentModel$formula, currentModel$model)
  nObs <- max(currentModel$data$obj)

  matListSumObs <- vector("list", nObs)
  for( u in 1:nObs ){
    
    n_i <- max(currentModel$data[currentModel$data$obj==u, "timeInt"])
    matListSum <- vector("list", n_i)
    for( v in 1:n_i ){
      
      if(v==1){
        
        matListSum[[v]] <- hessianLogitLikGEE(k=k, i=u, l=v, modelEst=modelEst)
        
      } else{
        
        matListSum[[v]] <- matListSum[[v-1]] + hessianLogitLikGEE(k=k, i=u, l=v, modelEst=modelEst)
        
      }
      
    }
    
    if(u==1){
      
      matListSumObs[[u]] <- matListSum[[n_i]]
      
    } else{
      
      matListSumObs[[u]] <- matListSumObs[[u-1]] + matListSum[[n_i]]
      
    }
    
    
    
  }
  
  - matListSumObs[[nObs]] / nObs
  
}

#########################################################
# Common covariance of parameter estimates between events

#' GEE covariance of all events for discrete competing risks
#' 
#' @description Estimates covariance of estimated parameters of all competing events generalized estimation equation models using sandwich approach.
#' @param modelEst Discrete time competing risks GEE model prediction model ("class dCRGEE").
#' @return Returns symmetric matrix of rows and columns dimension "number of competing risks" * "number of regression parameters" ("numeric matrix").
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link{compRisksGEE}}, \code{\link{dataLongCompRisks}}, \code{\link{dataLongCompRisksTimeDep}}, 
#' \code{\link[geepack]{geeglm}}
#' @references 
#' \insertRef{minjungDiscComp}{discSurv}
#' @keywords survival
#' @examples
#' 
#' # Example with unemployment data
#' library(Ecdat)
#' data(UnempDur)
#' 
#' # Select subsample
#' SubUnempDur <- UnempDur [1:100, ]
#' 
#' # Estimate GEE models for all events
#' estGEE <- compRisksGEE(datShort = SubUnempDur, dataTransform = "dataLongCompRisks", 
#' corstr = "independence", formulaVariable =~ timeInt + age + ui + logwage * ui, 
#' eventColumns = c("censor1", "censor2", "censor3", "censor4"), timeColumn = "spell")
#' 
#' \dontrun{
#' # Estimate covariance matrix of estimated parameters and competing events
#' estCovar <- covarGEE(modelEst=estGEE)
#' estCovar
#' 
#' # Covariances of estimated parameters of one event equal the diagonal blocks
#' lengthParameters <- length(estGEE[[1]]$coefficients)
#' noCompEvents <- length(estGEE)
#' meanAbsError <- rep(NA, noCompEvents)
#' for( k in 1:noCompEvents ){
#'   
#'   relInd <- (1 + (k-1) * lengthParameters) : (k * lengthParameters)
#'   meanAbsError[k] <- mean(abs(estCovar[relInd, relInd] - estGEE[[k]]$geese$vbeta))
#'   
#' }
#' mean(meanAbsError) 
#' # -> Covariance estimates within each event are equal to diagonal blocks in 
#' # complete covariance matrix with very small differences due to numerical accuracy.
#' }
#' 
#' @export covarGEE
covarGEE <- function(modelEst){
  
  # 1. Matrix A
  nObs <- max(modelEst[[1]]$data$obj)
  k_max <- length(modelEst)
  pMax <- length(modelEst[[1]]$coefficients)
  matA <- matrix(0, nrow=k_max*pMax, ncol=k_max*pMax)
  for( k in 1:k_max ){
    
    tempMat <- H_k(k=k, modelEst=modelEst)
    relInd <- ( 1 + (k-1) * pMax):(k * pMax)
    matA[relInd, relInd] <- solve(tempMat)
    
  }
  
  # 2. Matrix B
  matBsum <- matrix(0, nrow=k_max * pMax, ncol=k_max * pMax)
  for( i in 1:nObs ){
    
    uVec <- rep(NA, k_max * pMax)
    for( k in 1:k_max ){
      
      relInd <- ( 1 + (k-1) * pMax):(k * pMax)
      uVec[relInd] <- U_k_i(k=k, i=i, modelEst=modelEst)
      
    }
    
    matBsum <- matBsum + outer(uVec, uVec)
    
  }
  matB <- matBsum / nObs
  
  # 3. Output
  finalMat <- matA %*% matB %*% matA / nObs
  namesTemp <- c(sapply(1:k_max, function(x) paste(
    names(modelEst)[[x]], "_", names(modelEst[[x]]$coefficients), sep="")))
  dimnames(finalMat) <- list(namesTemp, namesTemp)
  return(finalMat)
  
}

###############################################################################

#' Estimates Cumulative Incidence Function for Discrete Time Competing Risks Models
#' 
#' Estimates the cumulative incidence function of a discrete time competing risks model 
#' given covariates P(T <= t, event = k | x). 
#' 
#' @param hazards Estimated discrete hazard rates of all events ("numeric matrix"). 
#' Each column represents one event. The first column is assumed to contain the censoring case 
#' and the discrete hazards should only vary over time in each row. 
#' @param eventFocus Column that represent the discrete hazards of the primary event ("integer vector").
#' @return Returns cumulative incidence function of the primary event. 
#' If argument \emph{nonparCI} is set to TRUE, then a list is returned: 
#' The first element includes the cumulative incidence function. 
#' The second list element contains the lower and the third list element the upper bound of
#' the pointwise confidence intervals.
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @details The covariates set is required to be constant across rows. 
#' @seealso \code{\link{compRisksGEE}}, \code{\link{dataLongCompRisks}}, 
#' \code{\link{dataLongCompRisksTimeDep}}, \code{\link[geepack]{geeglm}}, 
#' @references 
#' \insertRef{minjungDiscComp}{discSurv}
#' @keywords survival
#' @examples
#' 
#' # Example with unemployment data
#' library(Ecdat)
#' data(UnempDur)
#' 
#' # Select subsample
#' SubUnempDur <- UnempDur [1:100, ]
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
#' # Estimate cumulative incidence function
#' cumInzGEE <- estCumInz(preds, eventFocus = 2)
#' cumInzGEE
#' 
#' @export estCumInz
estCumInz <- function(hazards, eventFocus){
  
  # Estimate cumulative incidence function
  allCauseHazard <- rowSums(hazards)
  P_T_k <- hazards[, eventFocus] * c(1, estSurv(allCauseHazard)[1:(dim(hazards)[1] - 1)] )
  cumInz1 <- cumsum(P_T_k)
  names(cumInz1) <- paste("P(T<=", 1:length(cumInz1), ", event)", sep = "")
  return(cumInz1)
}

##########################
# gumbel

# Description:
# Specifies the link function of gumbel distribution function in the contex of general, linear models
# Can be used in function glm

# Output
# Should be used in context with glm in the base package

#' Gumbel Link Function
#' 
#' Constructs the link function with gumbel distribution in approriate format
#' for use in generalized, linear models.
#' 
#' Insert this function into a binary regression model
#' 
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' 
#' Matthias Schmid \email{matthias.schmid@@imbie.uni-bonn.de}
#' @seealso \code{\link{glm}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv}
#' @keywords survival
#' @examples
#' 
#' # Example with copenhagen stroke study
#' library(pec)
#' data(cost)
#' head(cost)
#' 
#' # Take subsample and convert time to months
#' costSub <- cost [1:50, ]
#' costSub$time <- ceiling(costSub$time/30)
#' costLong <- dataLong(dataShort = costSub, timeColumn = "time", eventColumn = "status",
#' timeAsFactor=TRUE)
#' gumbelModel <- glm(formula = y ~ timeInt + diabetes, data = costLong, 
#' family = binomial(link = gumbel()))
#' 
#' # Estimate hazard given prevStroke and no prevStroke
#' hazPrevStroke <- predict(gumbelModel, newdata=data.frame(timeInt = factor(1:143), 
#' diabetes = factor(rep("yes", 143), levels = c("no", "yes"))), type = "response")
#' hazWoPrevStroke <- predict(gumbelModel, newdata = data.frame(timeInt = factor(1:143), 
#' diabetes=factor(rep("no", 143), levels = c("no", "yes"))), type = "response")
#' 
#' # Estimate survival function
#' SurvPrevStroke <- cumprod(1 - hazPrevStroke)
#' SurvWoPrevStroke <- cumprod(1 - hazWoPrevStroke)
#' 
#' # Example graphics of survival curves with and without diabetes
#' plot(x = 1:143, y = SurvWoPrevStroke, type = "l", xlab = "Months", 
#' ylab = "S (t|x)", las = 1, lwd = 2, ylim = c(0,1))
#' lines(x = 1:143, y = SurvPrevStroke, col = "red", lwd = 2)
#' legend("topright", legend = c("Without diabetes", "Diabetes"), 
#' lty = 1, lwd =2, col = c("black", "red"))
#' 
#' @export gumbel
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



#' Estimated Survival Function
#' 
#' Estimates the survival function S(T = t|x) based on estimated hazard rates.
#' The hazard rates may or may not depend on covariates. The covariates have to
#' be equal across all estimated hazard rates. Therefore the given hazard rates
#' should only vary over time.
#' 
#' The argument \emph{haz} must be given for all intervals [a_0, a_1), [a_1,
#' a_2), ..., [a_{q-1}, a_q), [a_q, Inf).
#' 
#' @param haz Estimated hazard rates ("numeric vector")
#' @return Estimated probabilities of survival ("numeric vector")
#' @note It is assumed that all time points up to the last 
#' theoretical interval [a_q, Inf) are available. If not already present, 
#' these can be added manually.
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link{estMargProb}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv}
#' @keywords survival
#' @examples
#' 
#' 
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
#' Fit <- glm(formula = y ~ timeInt + age + logwage, data=UnempLong, family = binomial())
#' 
#' # Estimate discrete survival function given age, logwage of first person
#' hazard <- predict(Fit, newdata = subset(UnempLong, obj == 1), type = "response")
#' SurvivalFuncCondX <- estSurv(c(hazard, 1))
#' SurvivalFuncCondX
#' 
#' @export estSurv
estSurv <- function (haz) {
  
  # Input checks
  if(!is.vector(haz)) {stop(
    "Argument *haz* is not a vector! Please specify a numeric vector of estimated hazard rates.")}
  if(!all(-sqrt(.Machine$double.eps) <= haz & 
          haz - 1 <= sqrt(.Machine$double.eps))) {
    stop("Argument *haz* is not a vector of probabilities! Please specify a numeric vector of estimated hazard rates.")}
  
  erg <- cumprod(1 - haz)
  names(erg) <- paste("S(T=", 1:length(haz), ")", sep="")
  return(erg)
}

#########################
# estSurvCens

# Description
# Estimates the censoring survival function using life table estimator

#' Estimated Survival Function of Censoring Process
#' 
#' Estimates the marginal survival function G(T=t) of the censoring process based on 
#' a life table estimator. Compatible with single event and competing risks data.
#' 
#' 
#' @param dataShort Data in original short format ("class data.frame").
#' @param timeColumn Name of column with discrete time intervals ("character
#' vector").
#' @param eventColumns Names of the event columns of \code{dataShort}("character vector"). In the
#' competing risks case the event columns have to be in dummy encoding format
#' ("numeric vector").
#' @return Named vector of estimated survival function of the censoring process
#' for all time points except the last theoretical interval.
#' @note In the censoring survival function the last time interval [a_q, Inf)
#' has the probability of zero.
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link{estSurv}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv}
#' @keywords survival
#' @examples
#' 
#' 
#' # Load unemployment data
#' library(Ecdat)
#' data(UnempDur)
#' 
#' # Select subsample
#' subUnempDur <- UnempDur [1:100, ]
#' 
#' ######################
#' # Single event example
#' 
#' # Estimate censoring survival function G(t)
#' estG <- estSurvCens(dataShort = subUnempDur, timeColumn = "spell", 
#' eventColumns = "censor1")
#' estG
#' 
#' #########################
#' # Competing risks example
#' 
#' # Estimate censoring survival function G(t)
#' estG <- estSurvCens(dataShort = subUnempDur, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor2", "censor3", "censor4"))
#' estG
#' 
#' 
#' @export estSurvCens
estSurvCens <- function(dataShort, timeColumn, eventColumns) {

  # # Expand short data to long competing risks format
  # dataLongFormat  <- dataLongCompRisks (dataShort = dataShort, timeColumn = timeColumn, eventColumns = eventColumns, 
  #                    timeAsFactor = timeAsFactor)
  # 
  # # Construct censoring data
  # dataShortLongCens <- dataCensoring (dataShortLong = dataLongFormat, respColumn = "y", 
  #                                   timeColumn = "timeInt")

  # Convert to censoring format
  dataShortCens <- dataCensoring (dataShort = dataShort, eventColumns = eventColumns,
                                     timeColumn = timeColumn)

  # Estimate nonparametric survival function of censoring variable
  tempLifeTab <- lifeTable(dataShort = dataShortCens, timeColumn = "timeCens", 
                            eventColumn = "yCens")
  preG <- tempLifeTab [[1]] [, "S"]
  GT <- c(1, preG)
  names(GT) <- paste("t=", 0:length(preG), sep="")
  return(GT)
}

#########################
# estMargProb

# Description
# Estimates the discrete marginal probability P(T = t) given discrete hazard rates
# The hazard rates may or not depend on covariates, but the covariates have to be equal for all estimates of the hazard rate!
# Therefore the hazard rates only depend on time. 

# Input
# haz: Estimated hazard rates, which only depend on a time interval (numeric vector). 
# It is assumed that the estimates contain the intervals [0, a_1), [a_1, a_2), ..., [a_{q-1}, a_{q}) and the last interval is computed

# Output
# Survival function over all intervals [0, a_1), [a_1, a_2), ..., [a_{q-1}, a_{q}), [a_q, a_{\infty}



#' Estimated Marginal Probabilities
#' 
#' Estimates the marginal probability P(T=t|x) based on estimated discrete hazards.
#' The discrete hazards may or may not depend on covariates. The covariates have to
#' be equal across all estimated hazard rates. Therefore the given discrete hazards
#' should only vary over time.
#' 
#' The argument \emph{hazards} must be given for all intervals [a_0, a_1), [a_1,
#' a_2), ..., [a_{q-1}, a_q), [a_q, Inf).
#' 
#' @param hazards Estimated discrete hazards ("numeric vector")
#' @return Estimated marginal probabilities ("numeric vector")
#' @note It is assumed that all time points up to the last interval [a_q, Inf)
#' are available. If not already present, these can be added manually.
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link{estSurv}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv}
#' @keywords survival
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
#' UnempLong <- dataLong(dataShort = subUnempDur, timeColumn = "spell", eventColumn = "censor1")
#' head(UnempLong)
#' 
#' # Estimate binomial model with logit link
#' Fit <- glm(formula = y ~ timeInt + age + logwage, data = UnempLong, family = binomial())
#' 
#' # Estimate discrete survival function given age, logwage of first person
#' hazard <- predict(Fit, newdata = subset(UnempLong, obj == 1), type = "response")
#' 
#' # Estimate marginal probabilities given age, logwage of first person
#' MarginalProbCondX <- estMargProb (c(hazard, 1))
#' MarginalProbCondX
#' sum(MarginalProbCondX)==1 # TRUE: Marginal probabilities must sum to 1!
#' 
#' @export estMargProb
estMargProb <- function (hazards) {
  
  # Input checks
  if(!is.vector(hazards)) {stop("Argument *hazards* is not a vector! Please specify a numeric vector of estimated hazard rates.")}
  if(!all(-sqrt(.Machine$double.eps) <= hazards & 
          hazards - 1 <= sqrt(.Machine$double.eps))) {
    stop("Argument *hazards* is not between 0 and 1! Please specify a numeric vector of probabilities.")
  }
  if(length(hazards) == 0) {stop("Argument *hazards* is empty! Please specify a numeric vector of probabilities.")}
  
  if(length(hazards) > 1) {
    EstSurv <- estSurv(hazards)
    EstProb <- c(hazards[1], sapply(2:length(hazards), function (x) hazards [x] * EstSurv [x - 1]))
    names(EstProb) <- paste("P(T=", 1:length(hazards), ")", sep = "")
  }
  else {
    EstProb <- 1
    names(EstProb) <- paste("P(T=", 1, ")", sep = "")
  }

  return(EstProb)
}
