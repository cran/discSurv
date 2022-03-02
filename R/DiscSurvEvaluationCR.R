#' Discrete concordance index for competing risks 
#' 
#' Estimates the discrete concordance index in the case of competing risks.
#' 
#' @param markers Predictions on the test data with model fitted on training data ("numeric matrix").
#' Predictions are stored in the rows and the number of columns equal to the number of events. 
#' @param testTime New time intervals in the test data ("integer vector").
#' @param testEvents New event indicators (0 or 1) in the test data ("binary matrix"). Number of columns are
#' equal to the number of events.  
#' @param trainTime Time intervals in the training data ("integer vector").
#' @param trainEvents Event indicators (0 or 1) in the training data ("binary matrix"). Number of columns are
#' equal to the number of events.
#' @return Value of discrete concordance index between zero and one ("numeric vector").
#' @note It is assumed that all time points up to the last observed interval
#' [a_{q-1}, a_q) are available.
#' @author  Moritz Berger \email{moritz.berger@@imbie.uni-bonn.de} \cr \url{https://www.imbie.uni-bonn.de/personen/dr-moritz-berger/}
#' @seealso \code{\link{cIndex}}
#' @references 
#' \insertRef{heyardValCompRisks}{discSurv}
#' @keywords discrimination discrete_survival competing_risks
#' @examples
#' 
#' ##################################################
#' # Example with unemployment data and prior fitting
#' 
#' library(Ecdat)
#' data(UnempDur)
#' summary(UnempDur$spell)
#' # Extract subset of data
#' set.seed(635)
#' IDsample <- sample(1:dim(UnempDur)[1], 100)
#' UnempDurSubset <- UnempDur [IDsample, ]
#' set.seed(-570)
#' TrainingSample <- sample(1:100, 75)
#' UnempDurSubsetTrain <- UnempDurSubset [TrainingSample, ]
#' UnempDurSubsetTest <- UnempDurSubset [-TrainingSample, ]
#' 
#' # Convert to long format
#' UnempDurSubsetTrainLong <- dataLongCompRisks(dataShort = UnempDurSubsetTrain, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor4"), timeAsFactor = TRUE)
#' 
#' # Estimate continuation ratio model with logit link
#' vglmFit <- VGAM::vglm(formula = cbind(e0, e1, e2) ~ timeInt + age + logwage, 
#' data = UnempDurSubsetTrainLong, family=VGAM::multinomial(refLevel = "e0"))
#' 
#' gamFitPreds <- VGAM::predictvglm(vglmFit , newdata = cbind(UnempDurSubsetTest, 
#' timeInt = as.factor(UnempDurSubsetTest$spell)))
#' 
#' # Evaluate C-Index based on short data format
#' cIndexCompRisks(markers = gamFitPreds, 
#' testTime = UnempDurSubsetTest$spell, 
#' testEvents = UnempDurSubsetTest[, c("censor1", "censor4")], 
#' trainTime = UnempDurSubsetTrain$spell, 
#' trainEvents = UnempDurSubsetTrain[, c("censor1", "censor4")])
#' 
#' @export cIndexCompRisks
cIndexCompRisks <- function(markers, testTime, testEvents, trainTime, trainEvents){
  Cindices <- numeric(ncol(markers))
  for(i in 1:ncol(markers)){
    Cindices[i] <- cIndex(markers[,i], testTime, testEvents[, i], trainTime, trainEvents[, i])
  }
  weights <- colSums(testEvents) / sum(colSums(testEvents))
  Cindex  <- sum(Cindices*weights)
  return(Cindex)
} 

#################################################
#' Prediction Error Curves for Competing Risks 
#' 
#' Estimates prediction error curves for discrete survival competing risks models
#' 
#' @param testPreds Predictions on the test data with model fitted on training data ("numeric matrix").
#' Predictions are stored in the rows and the number of columns equal the number of events.
#' @param testDataShort Test data in short format ("class data.frame").
#' @param trainDataShort Train data in short format ("class data.frame").
#' @param timeColumn Character giving the column name of the observed times("character vector").
#' @param eventColumns Character vector giving the column names of the event indicators (excluding censoring column) ("character vector").
#' @param tmax Gives the maximum time interval for which prediction errors are
#' calculated ("integer vector"). It must not be higher than the maximum observed time in the
#' training data.
#' @return Calculated prediction errors for each competing event. Array with one matrix per competing event, 
#' with the predictions in the rows and the time points in the columns. 
#' @author Moritz Berger \email{moritz.berger@@imbie.uni-bonn.de} \cr \url{https://www.imbie.uni-bonn.de/personen/dr-moritz-berger/}
#' @seealso \code{\link{intPredErrCompRisks}}, \code{\link{predErrCurve}} 
#' @references 
#' \insertRef{heyardValCompRisks}{discSurv}
#' @keywords prediction_error competing_risks discrete_survival
#' @examples
#' 
#'###########################
#'# Example unemployment data
#'library(Ecdat)
#'data(UnempDur)
#'
#'# Select subsample
#'selectInd1 <- 1:200
#'selectInd2 <- 201:400
#'trainSet <- UnempDur[which(UnempDur$spell %in% (1:10))[selectInd1], ]
#'testSet <- UnempDur[which(UnempDur$spell %in% (1:10))[selectInd2], ]
#'
#'# Convert to long format 
#'trainSet_long <- dataLongCompRisks(dataShort=trainSet, timeColumn="spell", 
#'eventColumns=c("censor1", "censor4"), timeAsFactor=TRUE)
#'tmax          <- max(trainSet$spell)
#'testSet_long <- dataLongCompRisks(dataShort=testSet, timeColumn="spell", 
#'eventColumns=c("censor1", "censor4"), aggTimeFormat = TRUE, lastTheoInt=tmax,
#'timeAsFactor=TRUE)
#'
#'# Estimate continuation ratio model with logit link
#'vglmFit <- VGAM::vglm(formula=cbind(e0, e1, e2) ~ timeInt + age + logwage, 
#'data=trainSet_long, family=VGAM::multinomial(refLevel="e0"))
#'
#'# Calculate predicted hazards
#'predHazards <- VGAM::predictvglm(vglmFit, newdata=testSet_long, type="response")
#'
#'# Compute prediction error 
#'predErrCompRisks(testPreds=predHazards[,-1], testSet, trainSet, "spell", 
#'c("censor1", "censor4"), tmax)
#'
#' 
#' @export predErrCompRisks
predErrCompRisks <- function(testPreds, 
                             testDataShort, 
                             trainDataShort,
                             timeColumn, 
                             eventColumns, 
                             tmax=NULL){
  
  if(is.null(tmax)){
    tmax <- max(testDataShort[,timeColumn])
  }
  if(tmax > max(trainDataShort[,timeColumn])){
    stop("Argument *tmax* is higher than the latest observed interval in training data.")
  }
  
  # censoring distribution
  formC    <- formula(paste("y~","timeInt"))
  datC     <- dataCensoring(dataShort=trainDataShort, timeColumn=timeColumn, eventColumns=eventColumns)
  if(max(trainDataShort[,timeColumn]) > max(datC[,"timeCens"])){
    stop("No censored observations in all observed intervals. The censoring distribution can not be estimated for all intervals in training data.")
  }
  datLongC <- dataLong(dataShort=datC, timeColumn="timeCens", eventColumn="yCens", timeAsFactor = TRUE)
  modC     <- glm(formC, data=datLongC, family=binomial(link="logit"))
  predC    <- predict(modC, type="response", newdata=data.frame("timeInt"=factor(1:max(datC[,"timeCens"]))))
  Ghat     <- c(1, cumprod(1-predC))
  ipw      <- t(sapply(1:nrow(testDataShort), function(j) any(testDataShort[j,eventColumns]>0)*1*((testDataShort[j,timeColumn]<=c(1:max(trainDataShort[,timeColumn])))*1)/Ghat[testDataShort[j,timeColumn]]+
                         (testDataShort[j,timeColumn]>(1:max(trainDataShort[,timeColumn])))*1/Ghat[-1]))
  
  l_hats <- array(NA, dim=c(nrow(testDataShort), tmax, ncol(testPreds)))
  for(i in 1:ncol(testPreds)){
    l_hats[,,i] <- matrix(testPreds[,i], nrow=nrow(testDataShort), ncol=tmax, byrow=T)
  }
  l_hat  <- apply(l_hats, c(1,2), sum)
  S_hat  <- t(apply(1-l_hat, 1, cumprod))
  F_hat  <- array(NA, dim=c(nrow(testDataShort), tmax, ncol(testPreds)))
  for(i in 1:ncol(testPreds)){
    F_hat[,,i] <- t(sapply(1:nrow(testDataShort), function(j) cumsum(c(1,S_hat[j,1:(tmax-1)])*l_hats[j,,i])))
  }
  PE_hat <- array(NA, dim=c(nrow(testDataShort), tmax, ncol(testPreds)), dimnames=list(1:nrow(testDataShort), 1:tmax, eventColumns))
  for(i in 1:ncol(testPreds)){
    PE_hat[,,i] <- ipw[,1:tmax]*t(sapply(1:nrow(testDataShort), function(j) (c(testDataShort[j,timeColumn]<=c(1:tmax) & testDataShort[j,eventColumns[i]]==1)*1- F_hat[j,,i])^2))
  }
  return(PE_hat)
}

#################################################
#' Integrated Prediction Error for Competing Risks 
#' 
#' Estimates integrated prediction errors of arbitrary prediction models in the case of competing risks. 
#' 
#' @param testPreds Predictions on the test data with model fitted on training data. Must be a matrix, 
#' with the predictions in the rows and the number of columns equal to the number of events.
#' @param testDataShort Test data in short format.
#' @param trainDataShort Train data in short format.
#' @param timeColumn Character giving the column name of the observed times.
#' @param eventColumns Character vector giving the column names of the event indicators (excluding censoring column).
#' @param tmax Gives the maximum time interval for which prediction errors are
#' calculated. It must not be higher than the maximum observed time in the
#' training data.
#' @return Integrated prediction errors for each competing event. Matrix, with the integrated predictions in the rows 
#' and the number of columns equal to the number of events. 
#' @author Moritz Berger \email{moritz.berger@@imbie.uni-bonn.de} \cr \url{https://www.imbie.uni-bonn.de/personen/dr-moritz-berger/}
#' @seealso \code{\link{predErrCompRisks}}, \code{\link{predErrCurve}} 
#' @references 
#' \insertRef{heyardValCompRisks}{discSurv}
#' @keywords competing_risks discrete_survival
#' @examples
#' 
#'###########################
#'# Example unemployment data
#'library(Ecdat)
#'data(UnempDur)
#'
#'# Select subsample
#'selectInd1 <- 1:200
#'selectInd2 <- 201:400
#'trainSet <- UnempDur[which(UnempDur$spell %in% (1:10))[selectInd1], ]
#'testSet <- UnempDur[which(UnempDur$spell %in% (1:10))[selectInd2], ]
#'
#'# Convert to long format 
#'trainSet_long <- dataLongCompRisks(dataShort=trainSet, timeColumn="spell", 
#'eventColumns=c("censor1", "censor4"), timeAsFactor=TRUE)
#'tmax          <- max(trainSet$spell)
#'testSet_long <- dataLongCompRisks(dataShort=testSet, timeColumn="spell", 
#'eventColumns=c("censor1", "censor4"), aggTimeFormat = TRUE, lastTheoInt=tmax,
#'timeAsFactor=TRUE)
#'
#'# Estimate continuation ratio model with logit link
#'vglmFit <- VGAM::vglm(formula=cbind(e0, e1, e2) ~ timeInt + age + logwage, 
#'data=trainSet_long, family=VGAM::multinomial(refLevel="e0"))
#'
#'# Calculate predicted hazards
#'predHazards <- VGAM::predictvglm(vglmFit, newdata=testSet_long, type="response")
#'
#'# Compute integrated prediction error 
#'intPredErrCompRisks(testPreds=predHazards[,-1], testSet, trainSet, "spell", 
#'c("censor1", "censor4"), tmax)
#'
#'
#' @importFrom VGAM vglm 
#' @export intPredErrCompRisks
intPredErrCompRisks <- function(testPreds, 
                                testDataShort,
                                trainDataShort,
                                timeColumn, 
                                eventColumns,
                                tmax=NULL){
  
  if(is.null(tmax)){
    tmax <- max(testDataShort[,timeColumn])
  }
  if(tmax > max(trainDataShort[,timeColumn])){
    stop("Argument *tmax* is higher than the latest observed interval in training data.")
  }
  
  # marginal probabilities 
  trainDataLong <- dataLongCompRisks(dataShort=trainDataShort, timeColumn=timeColumn, 
                                     eventColumns=eventColumns, timeAsFactor = TRUE)
  form0 <- paste0("cbind(", paste0("e", 0:ncol(testPreds), collapse=","),")~timeInt")
  form0 <- as.formula(form0)
  margFit <- VGAM::vglm(form0, data=trainDataLong, family=VGAM::multinomial(refLevel="e0"))
  #predMargData <- data.frame(timeInt=factor(1:max(trainDataShort[,timeColumn])))
  
  predMargData <- dataLongCompRisks(dataShort=testDataShort, timeColumn=timeColumn, 
                                    eventColumns=eventColumns, aggTimeFormat = TRUE, lastTheoInt=tmax,
                                    timeAsFactor=TRUE)
  
  margHaz <- VGAM::predictvglm(margFit, predMargData, type="response")[,-1]
  overallHaz <- rowSums(margHaz)
  Smarg   <- cumprod(1-overallHaz)
  margProbs <- sapply(1:ncol(testPreds), function(j) margHaz[,j]*c(1,Smarg[1:(max(trainDataShort[,timeColumn])-1)]))
  
  # predition error per time point 
  predErr <- predErrCompRisks(testPreds, trainDataShort, testDataShort, timeColumn, eventColumns, tmax)
  
  # integrated error 
  intPredErr <- matrix(NA, nrow=nrow(testDataShort), ncol=ncol(testPreds))
  colnames(intPredErr) <- eventColumns
  for(i in 1:ncol(testPreds)){
    intPredErr[,i] <- sapply(1:nrow(testDataShort), function(j) sum(predErr[j,,i]*margProbs[1:tmax,i]))
  }
  return(intPredErr)
}
