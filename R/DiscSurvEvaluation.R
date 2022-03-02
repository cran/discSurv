########################
#' Calibration Plots 

#'@description Calibration plot based on predictions. Overall root mean squared error (RMSE) of 
#'predicted and observed discrete hazards is calculated. 
#' 
#'@param testPreds Predictions on the validation data with model fitted on training data ("numeric vector").
#'@param testDataLong Validation data set in long format ("class data.frame").
#'@param weights optional vector of weights ("numeric vector"). The length of weights must be equal to the number of observations 
#'of the validation data set. 
#'@param K Number of subsets for plotting ("integer vector").
#'@param event Column names of the event to be considered for plotting (only in case of cause-specific hazards) ("character vector").
#'@param ... Additional arguments passed to \code{\link{plot}}.
#'
#'@return Calibration plot 
#'@author Moritz Berger \email{moritz.berger@@imbie.uni-bonn.de} \cr \url{https://www.imbie.uni-bonn.de/personen/dr-moritz-berger/}
#'@seealso \code{\link{estRecal}}, \code{\link{dataLong}}, \code{\link{dataLongCompRisks}}, \code{\link{dataLongSubDist}}
#'@references 
#'\insertRef{bergerTutorial}{discSurv} \cr\cr \insertRef{heyardValCompRisks}{discSurv} \cr\cr \insertRef{bergerAssessing}{discSurv}
#'@keywords validation discrete_survival
#'@examples
#' ####################
#' # Data preprocessing
#' 
#' # Example unemployment data
#' library(Ecdat)
#' data(UnempDur)
#' 
#' # Select subsample
#' selectInd1 <- 1:100
#' selectInd2 <- 101:200
#' trainSet <- UnempDur[which(UnempDur$spell %in% (1:10))[selectInd1], ]
#' valSet <- UnempDur[which(UnempDur$spell %in% (1:10))[selectInd2], ]
#' 
#' ####################
#' # One event
#'
#' # Convert to long format
#' trainSet_long <- dataLong(dataShort = trainSet, timeColumn = "spell", eventColumn = "censor1")
#' valSet_long <- dataLong(dataShort = valSet, timeColumn = "spell", eventColumn = "censor1")
#' 
#' # Estimate continuation ratio model with logit link
#' glmFit <- glm(formula = y ~ timeInt + age + logwage, data = trainSet_long, family = binomial())
#' 
#' # Calculate predicted hazards
#' predHazards <- predict(glmFit, newdata = valSet_long, type = "response")
#' 
#' # Calibration plot
#' calibrationPlot(predHazards, testDataLong = valSet_long)
#' 
#' ############################
#' # Two cause specific hazards 
#' 
#' # Convert to long format
#' trainSet_long <- dataLongCompRisks(dataShort = trainSet, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor4"))
#' valSet_long <- dataLongCompRisks(dataShort = valSet, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor4"))
#' 
#' # Estimate continuation ratio model with logit link
#' vglmFit <- VGAM::vglm(formula = cbind(e0, e1, e2) ~ timeInt + age + logwage, data = trainSet_long, 
#' family = VGAM::multinomial(refLevel = "e0"))
#' 
#' # Calculate predicted hazards
#' predHazards <- VGAM::predictvglm(vglmFit, newdata = valSet_long, type = "response")
#' 
#' # Calibration plots
#' calibrationPlot(predHazards, testDataLong = valSet_long)
#' calibrationPlot(predHazards, testDataLong = valSet_long, event = "e2")
#' 
#' ###############################
#' # Subdistribution hazards model
#' 
#' # Convert to long format
#' trainSet_long <- dataLongSubDist(dataShort = trainSet, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor4"), eventFocus = "censor1")
#' valSet_long <- dataLongSubDist(dataShort = valSet, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor4"), eventFocus = "censor1")
#' 
#' # Estimate continuation ratio model with logit link
#' glmFit <- glm(formula = y ~ timeInt + age + logwage, data = trainSet_long, 
#' family = binomial(), weights = trainSet_long$subDistWeights)
#' 
#' # Calculate predicted  hazards 
#' predHazards <- predict(glmFit, newdata = valSet_long, type = "response")
#' 
#' # Calibration plot 
#' calibrationPlot(predHazards, testDataLong = valSet_long, weights = valSet_long$subDistWeights)
#' 
#' @export calibrationPlot
calibrationPlot <- function(testPreds, testDataLong, weights = NULL, 
                            K = 10, 
                            event = "e1", ...){
  
  
  if( is.null(dim(testPreds)) ){
    
    quantile_hat <- quantile(testPreds, probs = seq(0, 1, length.out = K + 1))
    subset_ID    <- cut(testPreds, breaks = quantile_hat)
    subsets_y    <- split(testDataLong$y, subset_ID)
    subsets_l    <- split(testPreds, subset_ID)
    if(is.null(weights)){
      observedProbs    <- sapply(subsets_y, mean)
      predictedHazards <- sapply(subsets_l, mean)
    } else{
      subsets_w    <- split(weights, subset_ID)
      observedProbs  <- sapply(1:K, function(j) (sum((subsets_y[[j]] * subsets_w[[j]])) + 1)/(sum(subsets_w[[j]]) + 2))
      predictedHazards  <- sapply(1:K, function(j) (sum((subsets_l[[j]] * subsets_w[[j]])) + 1)/(sum(subsets_w[[j]]) + 2))
    }
    xmax <- ymax <- max(c(observedProbs, predictedHazards))
    xlim <- c(0, xmax)
    ylim <- c(0, ymax)
    plot(x = predictedHazards, y = observedProbs, pch = 19, xlim = xlim, ylim = ylim, ...)
    lines(x = c(0, 1), y = c(0, 1), lty = "dashed")
    
  } else{
    
    lambda_hat   <- testPreds[,event]
    quantile_hat <- quantile(lambda_hat, probs = seq(0, 1, length.out = K + 1))
    subset_ID    <- cut(lambda_hat, breaks = quantile_hat)
    subsets_y    <- split(testDataLong[,event], subset_ID)
    observedProbs  <- sapply(subsets_y, mean)
    subsets_l    <- split(lambda_hat, subset_ID)
    predictedHazards  <- sapply(subsets_l, mean)
    xmax <- ymax <- max(c(observedProbs, predictedHazards))
    xlim <- c(0, xmax)
    ylim <- c(0, ymax)
    plot(x = predictedHazards, y = observedProbs, pch = 19, xlim = xlim, ylim = ylim, ...)
    lines(x = c(0, 1), y = c(0, 1), lty = "dashed")
    
  }
  
  RMSE <- sqrt(mean((predictedHazards - observedProbs)^2))
  names(RMSE) <- "RMSE"
  return(RMSE)
  
}

#####################################
#' Calibration Plots for Model Object 

#'@description Calibration plot based on discrete hazard, discrete subdistribution hazard
#' or discrete cause-specific hazards model. 
#' 
#'@param predModel A model object of class \code{"glm"} or \code{"vglm"} ("class glm,vglm").
#'@param type One out of "Hazard", "SubDistHazard", "CauseSpecHazards" ("character vector").
#'@param newdataLong Optional validation data set in long format ("class data.frame").
#'@param K Number of subsets for plotting("integer vector").
#'@param event Character giving the column names of the event to be considered for plotting (only for type="CauseSpecHazards")("character vector").
#'@param ... additional arguments passed to \code{\link{plot}}.
#'
#'@return Calibration plot 
#'@author @author Moritz Berger <moritz.berger@imbie.uni-bonn.de> \cr \url{http://www.imbie.uni-bonn.de/personen/dr-moritz-berger/}
#'@seealso \code{\link{estRecal}}, \code{\link{dataLong}}, \code{\link{dataLongCompRisks}}, \code{\link{dataLongSubDist}}
#'@references 
#'\insertRef{bergerTutorial}{discSurv} \cr\cr \insertRef{heyardValCompRisks}{discSurv} \cr\cr \insertRef{bergerAssessing}{discSurv}
#'@keywords internal
#'@examples
#' ####################
#' # Data preprocessing
#' 
#' # Example unemployment data
#' library(Ecdat)
#' data(UnempDur)
#' 
#' # Select subsample
#' selectInd1 <- 1:100
#' selectInd2 <- 101:200
#' trainSet <- UnempDur[which(UnempDur$spell %in% (1:10))[selectInd1], ]
#' valSet <- UnempDur[which(UnempDur$spell %in% (1:10))[selectInd2], ]
#' 
#' ####################
#' # One event
#'
#' # Convert to long format
#' trainSet_long <- dataLong(dataShort = trainSet, timeColumn = "spell", eventColumn = "censor1")
#' valSet_long <- dataLong(dataShort = valSet, timeColumn = "spell", eventColumn = "censor1")
#' 
#' # Estimate continuation ratio model with logit link
#' glmFit <- glm(formula = y ~ timeInt + age + logwage, data = trainSet_long, family = binomial())
#' 
#' # Calibration plot
#' calibrationPlotModel(glmFit, newdataLong = valSet_long)
#' 
#' ############################
#' # Two cause specific hazards 
#' 
#' # Convert to long format
#' trainSet_long <- dataLongCompRisks(dataShort = trainSet, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor4"))
#' valSet_long <- dataLongCompRisks(dataShort = valSet, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor4"))
#' 
#' # Estimate continuation ratio model with logit link
#' vglmFit <- VGAM::vglm(formula = cbind(e0, e1, e2) ~ timeInt + age + logwage, data = trainSet_long, 
#' family = VGAM::multinomial(refLevel = "e0"))
#' 
#' # Calibration plots
#' calibrationPlotModel(vglmFit, newdataLong = valSet_long, type = "CauseSpecHazard")
#' calibrationPlotModel(vglmFit, newdataLong = valSet_long, type = "CauseSpecHazard", event = "e2")
#' 
#' ###############################
#' # Subdistribution hazards model
#' 
#' # Convert to long format
#' trainSet_long <- dataLongSubDist(dataShort = trainSet, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor4"), eventFocus = "censor1")
#' valSet_long <- dataLongSubDist(dataShort = valSet, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor4"), eventFocus = "censor1")
#' 
#' # Estimate continuation ratio model with logit link
#' glmFit <- glm(formula = y ~ timeInt + age + logwage, data = trainSet_long, 
#' family = binomial(), weights = trainSet_long$subDistWeights)
#' 
#' # Calibration plot 
#' calibrationPlotModel(glmFit, newdataLong = valSet_long, type = "SubDistHazard")
#' 
#' @noRd
calibrationPlotModel <- function(predModel, 
                                 type = c("Hazard", "SubDistHazard", "CauseSpecHazards"), 
                                 newdataLong = NULL, 
                                 K = 10, 
                                 event = "e1", ...){
  
  type <- match.arg(type)
  
  if(type == "Hazard"){
    lambda_hat   <- predict(predModel, type = "response", newdata = newdataLong)
    quantile_hat <- quantile(lambda_hat, probs = seq(0, 1, length.out = K + 1))
    subset_ID    <- cut(lambda_hat, breaks = quantile_hat)
    if(is.null(newdataLong)){
      subsets    <- split(predModel$y, subset_ID)
    } else{
      subsets    <- split(newdataLong$y, subset_ID)
    }
    observedProbs  <- sapply(subsets, mean)
    subsets_l    <- split(lambda_hat, subset_ID)
    predictedHazards  <- sapply(subsets_l, mean)
    xmax <- ymax <- max(c(observedProbs, predictedHazards))
    xlim <- c(0, xmax)
    ylim <- c(0, ymax)
    plot(x = predictedHazards, y = observedProbs, pch = 19, xlim = xlim, ylim = ylim, ...)
    lines(x = c(0, 1), y = c(0, 1), lty = "dashed")
  }
  
  if(type == "CauseSpecHazards"){
    lambda_hat   <- VGAM::predict(predModel, type = "response", newdata = newdataLong)[, event]
    quantile_hat <- quantile(lambda_hat, probs = seq(0, 1, length.out = K + 1))
    subset_ID    <- cut(lambda_hat, breaks = quantile_hat)
    if(is.null(newdataLong)){
      subsets    <- split(predModel@y[,event], subset_ID)
    } else{
      subsets    <- split(newdataLong[,event], subset_ID)
    }
    observedProbs  <- sapply(subsets, mean)
    subsets_l    <- split(lambda_hat, subset_ID)
    predictedHazards  <- sapply(subsets_l, mean)
    xmax <- ymax <- max(c(observedProbs, predictedHazards))
    xlim <- c(0, xmax)
    ylim <- c(0, ymax)
    plot(x = predictedHazards, y = observedProbs, pch = 19, xlim = xlim, ylim = ylim, ...)
    lines(x = c(0, 1), y = c(0, 1), lty = "dashed")
  }
  
  if(type=="SubDistHazard"){
    lambda_hat   <- predict(predModel, type = "response", newdata = newdataLong)
    quantile_hat <- quantile(lambda_hat, probs = seq(0, 1, length.out = K + 1))
    subset_ID    <- cut(lambda_hat, breaks = quantile_hat)
    if(is.null(newdataLong)){
      subsets_y    <- split(predModel$y, subset_ID)
      subsets_w    <- split(predModel$weights, subset_ID)
    } else{
      subsets_y    <- split(newdataLong$y, subset_ID)
      subsets_w    <- split(newdataLong$subDistWeights, subset_ID)
    }
    observedProbs  <- sapply(1:K, function(j) (sum((subsets_y[[j]] * subsets_w[[j]])) + 1) / (sum(subsets_w[[j]]) + 2))
    subsets_l    <- split(lambda_hat, subset_ID)
    predictedHazards  <- sapply(1:K, function(j) (sum((subsets_l[[j]] * subsets_w[[j]])) + 1)/(sum(subsets_w[[j]]) + 2))
    xmax <- ymax <- max(c(observedProbs, predictedHazards))
    xlim <- c(0, xmax)
    ylim <- c(0, ymax)
    plot(x = predictedHazards, y = observedProbs, pch = 19, xlim = xlim, ylim = ylim, ...)
    lines(x = c(0, 1), y = c(0, 1), lty = "dashed")
  }
}

########################
# Logistic recalibration

# Cases
# One event
# No dependency time 
# Time dependent covariates

# Multiple events
# Competing risks
# Competing risks time dependency
# Subdistribution

#' Logistic recalibration based on linear predictors
#' 
#' Fits a logistic recalibration model to independent test data. It updates the intercept and slope parameters. 
#' It is assumed that the factor levels of time are equal in both training and validation data. 
#' Time dependent covariates, discrete cause specific competing risks and subdistribution hazards are also supported. 
#' 
#' @param testLinPred Calculated linear predictor on the validation data with model fitted on training data ("numeric vector").
#' @param testDataLong Validation data set in long format ("class data.frame"). 
#' @param weights Weights used in estimation of the logistic recalibration model ("numeric vector"). 
#' Default is no weighting (NULL).  
#' @return Continuation ratio model that calibrates estimated discrete hazards to new validation environment ("class glm, lm").
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @details Updates estimated hazards of discrete survival models to better adapt to a different environment. 
#' If there are substantial environment changes the predicted probabilities will differ between two environments.
#' Logistic recalibration may be used to improve the calibration of predicted probabilities by 
#' incorperating information from the existing model and data from the environment. This approach works 
#' for any survival prediction model with one event that provides linear predictors.
#' @seealso (Calibration plots links) \code{\link{dataLong}}, \code{\link{dataLongTimeDep}}, \code{\link{dataLongCompRisks}}, 
#' \code{\link{dataLongCompRisksTimeDep}}, \code{\link{dataLongSubDist}}
#' @references 
#' \insertRef{heyardValCompRisks}{discSurv}
#' @keywords survival
#' @examples
#' 
#' ####################
#' # Data preprocessing
#' 
#' # Example unemployment data
#' library(Ecdat)
#' data(UnempDur)
#' 
#' # Select subsample
#' selectInd1 <- 1:100
#' selectInd2 <- 101:200
#' trainSet <- UnempDur[which(UnempDur$spell %in% (1:10))[selectInd1], ]
#' valSet <- UnempDur[which(UnempDur$spell %in% (1:10))[selectInd2], ]
#' 
#' ####################
#' # One event
#'
#' # Convert to long format
#' trainSet_long <- dataLong(dataShort = trainSet, timeColumn = "spell", eventColumn = "censor1")
#' valSet_long <- dataLong(dataShort = valSet, timeColumn = "spell", eventColumn = "censor1")
#' 
#' # Estimate continuation ratio model with logit link
#' glmFit <- glm(formula = y ~ timeInt + age + logwage, data = trainSet_long, family = binomial())
#' 
#' # Calculate linear predictors on validation set
#' linPred <- predict(glmFit, newdata = valSet_long, type = "link")
#' 
#' # Estimate logistic recalibration model
#' recalModel <- estRecal(testLinPred = linPred, testDataLong = valSet_long)
#' summary(recalModel)
#' 
#' # Calibration plots
#' hazOrg <- predict(glmFit, newdata = valSet_long, type = "response")
#' hazRecal <- predict(recalModel, newdata = data.frame(linPred), type = "response")
#' 
#' # Before logistic recalibration
#' calibrationPlot(hazOrg, testDataLong = valSet_long)
#' # After logistic recalibration
#' calibrationPlot(hazRecal, testDataLong = valSet_long)
#' 
#' ############################
#' # Two cause specific hazards 
#' library(VGAM)
#' 
#' # Convert to long format
#' trainSet_long <- dataLongCompRisks(dataShort = trainSet, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor4"))
#' valSet_long <- dataLongCompRisks(dataShort = valSet, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor4"))
#' 
#' # Estimate continuation ratio model with logit link
#' vglmFit <- VGAM::vglm(formula = cbind(e0, e1, e2) ~ timeInt + age + logwage, data = trainSet_long, 
#' family = VGAM::multinomial(refLevel = "e0"))
#' 
#' # Calculate linear predictors on training and test set
#' linPred <- VGAM::predictvglm(vglmFit, newdata = valSet_long, type = "link")
#' 
#' # Estimate logistic recalibration model
#' recalModel <- estRecal(testLinPred = linPred, testDataLong = valSet_long)
#' recalModel
#' 
#' # Calibration plots
#' hazOrg <- predict(vglmFit, newdata = valSet_long, type = "response")
#' predDat <- as.data.frame(linPred)
#' names(predDat) <- recalModel@misc$colnames.x[-1]
#' hazRecal <- predictvglm(recalModel, newdata = predDat, type = "response")
#' 
#' # Before logistic recalibration
#' calibrationPlot(hazOrg, testDataLong = valSet_long, event = "e1")
#' calibrationPlot(hazOrg, testDataLong = valSet_long, event = "e2")
#' # After logistic recalibration
#' calibrationPlot(hazRecal, testDataLong = valSet_long, event = "e1")
#' calibrationPlot(hazRecal, testDataLong = valSet_long, event = "e2")
#' 
#' ###############################
#' # Subdistribution hazards model
#' 
#' # Convert to long format
#' trainSet_long <- dataLongSubDist(dataShort = trainSet, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor4"), eventFocus = "censor1")
#' valSet_long <- dataLongSubDist(dataShort = valSet, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor4"), eventFocus = "censor1")
#' 
#' # Estimate continuation ratio model with logit link
#' glmFit <- glm(formula = y ~ timeInt + age + logwage, data = trainSet_long, 
#' family = binomial(), weights = trainSet_long$subDistWeights)
#' 
#' # Calculate linear predictors on training and test set
#' linPred <- predict(glmFit, newdata = valSet_long, type = "link")
#' 
#' # Estimate logistic recalibration model
#' recalModel <- estRecal(testLinPred = linPred, testDataLong = valSet_long, 
#' weights = valSet_long$subDistWeights)
#' recalModel
#' 
#' # Calibration plots
#' hazOrg <- predict(glmFit, newdata = valSet_long, type = "response",
#' weights = valSet_long$subDistWeights)
#' hazRecal <- predict(recalModel, newdata = data.frame(linPred), type = "response",
#' weights = valSet_long$subDistWeights)
#' 
#' # Before logistic recalibration
#' calibrationPlot(hazOrg, testDataLong = valSet_long,
#' weights=valSet_long$subDistWeights)
#' # After logistic recalibration
#' calibrationPlot(hazRecal, testDataLong = valSet_long,
#' weights=valSet_long$subDistWeights)
#' 
#' @export estRecal
estRecal <- function(testLinPred, testDataLong, weights = NULL) {
  
  # Check dimension of linear predictor
  if( is.null(dim(testLinPred)) ){
    
    # Fit logistic re-calibration method with one event
    mergedDat <- data.frame(testDataLong, linPred = testLinPred)
    glmRecal <- glm(formula=y ~ linPred, data = mergedDat, family = binomial(), weights = weights)
    return(glmRecal)
    
  } else{
    
    # Fit logistic re-calibration method for competing risks
    mergedDat <- data.frame(testDataLong, testLinPred)
    names(mergedDat) <- gsub(".", "", names(mergedDat), fixed = TRUE)
    
    # Fit logistic recalibration model
    formulaCreate <- formula(paste("cbind(", paste(paste("e", 0:dim(testLinPred)[2], sep = ""), collapse = ","),
                                   ") ~ ", 
                                   paste(tail(names(mergedDat), dim(testLinPred)[2]), collapse = " + ", sep = "") ))
    vglmRecal <- vglm(formula = formulaCreate, data = mergedDat, 
                      family = multinomial(refLevel = "e0"), weights = weights)
    return(vglmRecal)
    
  }
}

#########################
# brierScore

# Description:
# Computes the Brier Score of person i, based on generalized, linear models given Survival and Censoring functions

# Steps
# 0. Extract training and test data
# 1. Convert training data to long format
# 2. Convert test data to long format
# 3. Convert response in training data to censoring variable
# 4. Fit survival model on training data in long format
# 5. Fit censoring model on training data in long format
# 6. Predict survival times on test data
# 7. Predict censoring times on test data
# 8. Estimate survival functions of survival times for each person (aggregate)
# 9. Estimate survival functions of censoring times for each person (aggregate)
# 10. Calculate brier score

# Input
# dataShort: Original data in short format. Must be of type data.frame
# trainIndices: List of integer Indices, which give the rows of *dataShort* as training data in cross validation
# survModelFormula: Formula of the survival process. First argument must be the univariate time (numeric) variable
# censModelFormula: Formula of the censoring process. First argument must be the univariate event (binary) variable 
# linkFunc: Gives the type of link used in the glm function. May be customized
# idColumn: Gives the column name of the identification numbers as character scalar. 
# Default NULL means, that each row equals one person (no repeated measurements)

# Output
# Numeric vector giving the brier scores of the persons available in all test subsets

#' Brier Score Calculation
#' 
#' Computes the brier score of each person based on a generalized, linear model
#' with cross validation.
#' 
#' At normal circumstances a covariate free model for the censoring weights is
#' sufficient.
#' 
#' @param dataShort Original data in short format ("class data.frame")
#' @param trainIndices List of indices used for training data sets ("class list").
#' E. g. if a 10-fold cross validation should be conducted, then the list
#' should have length 10 and in each element are the training indices ("integer vector").
#' @param survModelFormula Gives the specified relationship of discrete
#' response and covariates ("class formula"). The formula is designed, that the intercepts for
#' the time dependent base line hazards are always included. Therefore only
#' covariates should be given in this formula ("class "formula").
#' @param linkFunc Specifies the desired link function in use of generalized,
#' linear models ("character vector").
#' @param eventColumn Gives the column name of the event indicator (1=observed,
#' 0=censored) ("character vector").
#' @param idColumn Gives the column name of the identification number of each
#' person ("character vector").
#' Default NULL means, that each row equals one person (no repeated measurements).
#' @return Numeric vector of brier scores for each person in the original data.
#' @note It is assumed that all time points up to the last interval [a_q, Inf)
#' are available in short format. If not already present, these can be added
#' manually.
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' 
#' Matthias Schmid \email{matthias.schmid@@imbie.uni-bonn.de}
#' @seealso \code{\link{adjDevResidGlm}}, \code{\link{devResid}},
#' \code{\link{glm}}, \code{\link{predErrCurve}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv} \cr\cr
#' \insertRef{gerdsConsisEst}{discSurv}
#' @keywords internal
#' @examples
#' 
#' # Example with cross validation and unemployment data 
#' library(Ecdat)
#' library(caret)
#' data(UnempDur)
#' summary(UnempDur$spell)
#' 
#' # Extract subset of data
#' set.seed(635)
#' IDsample <- sample(1:dim(UnempDur)[1], 100)
#' UnempDurSubset <- UnempDur [IDsample, ]
#' head(UnempDurSubset)
#' range(UnempDurSubset$spell)
#' set.seed(7550)
#' CVfolds <- createFolds (y = UnempDurSubset$spell, returnTrain = TRUE, k = 2)
#' 
#' # Calculate brier score
#' tryBrierScore <- brierScore (dataShort = UnempDurSubset, trainIndices = CVfolds, 
#' survModelFormula = spell ~ age + logwage, linkFunc = "logit", 
#' eventColumn = "censor1", idColumn = NULL)
#' tryBrierScore
#' 
#' @noRd
brierScore <- function (dataShort, trainIndices, survModelFormula, linkFunc = "logit", eventColumn, idColumn = NULL) {
  
  # Input Checks
  if(!is.data.frame(dataShort)) {stop("Argument *dataShort* is not in the correct format! Please specify as data.frame object.")}
  if(!is.list(trainIndices)) {stop("Argument *trainIndices* is not in the correct format! Please specify a list.")}
  InputCheck1 <- all(sapply(1:length(trainIndices), function (x) is.integer(trainIndices [[x]])))
  if(!InputCheck1) {stop("Sublists of *trainIndices* are not all integer values! Please specify a list of integer Indices.")}
  if(!("formula" %in% class(survModelFormula))) {stop("*survModelFormula* is not of class formula! Please specify a valid formula, e. g. y ~ x + z.")}
  if(!any(names(dataShort) == eventColumn)) {stop("Argument *eventColumn* is not available in *dataShort*! Please specify the correct column name of the event indicator.")}
  
  # Help function
  B <- function(k) {
    probs <- estMargProb (LambdaSplit [[k]] [, "Lambda" ])
    # probs <- probs[-length(probs)]
    
    if(length(probs [-length(probs)]) != 0) {
      brierVec <- as.numeric(tail(LambdaSplit [[k]] [, eventColumn], 1) * (1 - tail(probs, 1))^2 + sum (probs [-length(probs)]))
    }
    else {
      brierVec <- as.numeric(tail(LambdaSplit [[k]] [, eventColumn], 1) * (1 - tail(probs, 1))^2)
    }
    return(brierVec)
  }
  
  RET <- vector("list", length(trainIndices))
  for(i in 1:length(trainIndices)) {
    
    # 0. Extract training, test data and responses
    TrainSet <- dataShort [trainIndices [[i]], ]
    if(length(trainIndices) != 1) {
      TestSet <- dataShort [-trainIndices [[i]], ]
    }
    else {
      TestSet <- TrainSet
    }
    
    # 1. Convert training data to long format
    if(!is.null(idColumn)) {
      TrainLong <- dataLongTimeDep(dataSemiLong = TrainSet, timeColumn = as.character(survModelFormula) [2], eventColumn = eventColumn, idColumn = idColumn)
    }
    else {
      TrainLong <- dataLong(dataShort = TrainSet, timeColumn = as.character(survModelFormula) [2], eventColumn = eventColumn)
    }
    
    # 2. Convert response in training data to censoring variable
    TrainLong <- dataCensoring(dataShort = TrainLong, timeColumn = "timeInt", shortFormat = FALSE)
    
    # 3. Convert test data to long format
    if(!is.null(idColumn)) {
      TestLong <- dataLongTimeDep(dataSemiLong = TestSet, timeColumn = as.character(survModelFormula) [2], eventColumn = eventColumn, idColumn = idColumn)
    }
    else {
      TestLong <- dataLong(dataShort = TestSet, timeColumn = as.character(survModelFormula) [2], eventColumn = eventColumn)
    }

    SurvnewFormula <- update(survModelFormula, y ~ timeInt + .)
    SurvFit <- glm (formula = SurvnewFormula, data = TrainLong, family = binomial(link = linkFunc), control = glm.control(maxit = 2500))

    # 6. Estimate survival curves on test data
    Check <- "error" %in% class(tryCatch(predict(SurvFit, TestLong, type = "response"), error = function (e) e))
    if(Check) {
      # Which columns are factors in test data?
      IndexFactor <- which(sapply(1:dim(TestLong)[2], function (x) is.factor(TestLong [, x])) == TRUE)
      # What are the levels of these factors?
      TestLevelsFactor <- sapply(IndexFactor, function (x) levels(TestLong [, x]))
      # First column does not count (censoring process)
      # What are the levels of the corresponding factors in the training data?
      TrainLevelsFactor <- sapply(IndexFactor, function (x) levels(TrainLong [, x + 1]))
      # Which levels of the test data exist in the training data factors?
      InLevelsFactor <- lapply(1:length(TestLevelsFactor), function (x) which((TestLevelsFactor [[x]] %in% TrainLevelsFactor [[x]]) == FALSE))
      ExcludeRows <- lapply (1:length(IndexFactor), function (j) which(TestLong [, IndexFactor [j]] %in% TestLevelsFactor [[j]] [InLevelsFactor [[j]]]))
      ExcludeRows <- do.call(c, ExcludeRows)
      TestLong <- TestLong [-ExcludeRows, ]
    }
    Lambda <- predict(SurvFit, TestLong, type = "response")
    LambdaSplit <- split(cbind(Lambda = Lambda, TestLong), TestLong$obj)
    RET [[i]] <- sapply(1:length(LambdaSplit), B)
  }

  # Output
  RET <- do.call(cbind, RET)
  RET <- rowMeans(RET)
  return(RET)
}

#############################
# tprUno
#############################

##############
# Description
# Estimates the predictive true positive rate (tpr) based on cross validation and generalized, linear models

#######
# Input
# timepoint: Discrete time interval given that the false positive rate is evaluated (integer scalar)
# dataShort: Original data. Should be in format data.frame()
# trainIndices: List of Indices from original data used for training (list of integer vectors). 
# The length of the list is equal to the number of cross valdiation samples
# survModelFormula: Formula of the survival model
# censModelFormula: Formula of the censoring model. Normally this is done without covariates
# linkFunc: Link function of the generalized, linear model see glm
# idColumn: Name of the column with identification numbers of persons. 
# Default NULL means, that each row equals one person (no repeated measurements).

# Output
# data.frame with columns
# cutoff: Cut off values of the linear predictor (numeric vector)
# fpr: False positive rate (numeric vector)



#' @name tprUno
#' @title True positive Rate Uno
#' 
#' @description Estimates the true positive rate based on Uno et al. to evaluate the
#' predictive accuracy of discrete generalized, linear survival models by cross
#' validation.
#' 
#' The formula \code{survModelFormula} has a specific structure: The
#' response on the left side of the formula is the time of the short data
#' format. On the right side are the covariates without time, e. g. Time ~ X1 +
#' X2 if there are only two covariates. The time will be added automatically.
#' 
#' 
#' @param timepoint Discrete time interval given that the false positive rate
#' is evaluated ("integer vector").
#' @param dataShort Original data in short format ("class data.frame").
#' @param trainIndices List of Indices from original data used for training
#' ("class list"). The length of the list is equal to the number of
#' cross validation samples.
#' @param survModelFormula Formula of the discrete survival model ("class formula"). It is used
#' in a generalized, linear model.
#' @param censModelFormula Formula of the censoring model ("class formula"). It is used in a
#' generalized, linear model. Usually this is done without covariates.
#' @param linkFunc Link function of the generalized, linear model ("character vector").
#' @param idColumn Name of the column with identification numbers of persons("character vector").
#' Default NULL means, that each row equals one person (no repeated
#' measurements).
#' @param timeAsFactor Should the time intervals be coded as factor ("logical vector")? Default is
#' to use factor. If the argument is false, the column is coded as numeric.
#' @return List with objects \itemize{ \item{Output} Data frame with two
#' columns: \emph{cutoff} gives the different marker values and \emph{tpr} the true
#' positive rates \item{Input} A list of given argument input values (saved for
#' reference). In addition there is the list element \code{orderMarker}, which
#' gives the indices of the marker values in increasing order. }
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' 
#' Matthias Schmid \email{matthias.schmid@@imbie.uni-bonn.de}
#' @seealso \code{\link[caret]{createDataPartition}}, \code{\link{glm}}
#' @references 
#' \insertRef{schmidDiscMeasure}{discSurv} \cr\cr
#' \insertRef{unoEvalPred}{discSurv} \cr\cr
#' \insertRef{heagertySurvROC}{discSurv}
#' @keywords internal
#' @examples
#' 
#' # Example with cross validation and unemployment data 
#' library(Ecdat)
#' library(caret)
#' data(UnempDur)
#' summary(UnempDur$spell)
#' 
#' # Extract subset of data
#' set.seed(635)
#' IDsample <- sample(1:dim(UnempDur)[1], 100)
#' UnempDurSubset <- UnempDur [IDsample, ]
#' head(UnempDurSubset)
#' range(UnempDurSubset$spell)
#' set.seed(7550)
#' CVfolds <- createFolds (y = UnempDurSubset$spell, returnTrain = TRUE, k = 2)
#' 
#' # Estimate true positive rate of time interval 7: 
#' # Correspondes to three and a half month duration (each interval is of length two weeks)
#' tryTPR <- tprUno (timepoint = 7, dataShort = UnempDurSubset, 
#' trainIndices = CVfolds, survModelFormula = spell ~ age + logwage, 
#' censModelFormula = censor1 ~ 1, linkFunc = "logit", idColumn = NULL)
#' tryTPR
#' plot(tryTPR)
#' 
#' @noRd
tprUno <- function(timepoint, dataShort, trainIndices, survModelFormula, censModelFormula, linkFunc = "logit", idColumn = NULL, timeAsFactor = TRUE) {
  
  # Input Checks
  if(length(timepoint) != 1 || !(timepoint == floor(timepoint))) {stop("Argument *timepoint* is not in the correct format! Please specify as integer scalar value.")}
  if(!is.data.frame(dataShort)) {stop("Argument *dataShort* is not in the correct format! Please specify as data.frame object.")}
  if(!is.list(trainIndices)) {stop("Argument *trainIndices* is not in the correct format! Please specify a list.")}
  InputCheck1 <- all(sapply(1:length(trainIndices), function (x) is.integer(trainIndices [[x]])))
  if(!InputCheck1) {stop("Sublists of *trainIndices* are not all integer values! Please specify a list of integer Indices.")}
  if(length(trainIndices) != 1) {
    InputCheck2 <- all(sort(as.numeric(do.call(c, lapply(trainIndices, function (x) setdiff(1:dim(dataShort) [1], x))))) == (1:dim(dataShort) [1]))
  }
  else {
    InputCheck2 <- all(trainIndices [[1]]==(1:dim(dataShort) [1]))
  }
  if(!InputCheck2) {stop("Argument *trainIndices* does not contain cross validation samples! Please ensure that the union of all test indices equals the indices of the complete data set.")}
  if(!("formula" %in% class(censModelFormula))) {stop("*censModelFormula* is not of class formula! Please specify a valid formula, e. g. yCens ~ 1")}
  if(!("formula" %in% class(survModelFormula))) {stop("*survModelFormula* is not of class formula! Please specify a valid formula, e. g. y ~ x")}
  if(!(any(names(dataShort) == idColumn) | is.null(idColumn))) {stop("Argument *idColumn* is not available in *dataShort*! Please specify the correct column name of the identification number.")}
  
  # Help function
  sens <- function(k) {
    
    sensNum <- sum((marker > k) * (newTime == timepoint) * newEvent / GT, na.rm = TRUE)
    sensDenom <- sum((newTime == timepoint) * newEvent / GT, na.rm = TRUE)
    
    if (sensDenom > 0)
      return(sensNum / sensDenom) else
        return(0)
  }
  
  # Loop across all training data sets
  markerList <- vector("list", length(trainIndices))
  ExcludeRowsCensList <- vector("list", length(trainIndices))
  ExcludeRowsdataShortList <- vector("list", length(trainIndices))
  oneMinuslambdaList <- vector("list", length(trainIndices))
  
  # Convert full sample to long format
  if(!is.null(idColumn)) {
    TrainLongFull <- dataLongTimeDep(dataSemiLong = dataShort, timeColumn = as.character(survModelFormula) [2], 
                                      eventColumn = as.character(censModelFormula) [2], 
                                      idColumn = idColumn, timeAsFactor = timeAsFactor)
  }
  else {
    TrainLongFull <- dataLong (dataShort = dataShort, timeColumn = as.character(survModelFormula) [2], 
                               eventColumn = as.character(censModelFormula) [2], timeAsFactor = timeAsFactor)
  }
  
  for(i in 1:length(trainIndices)) {
    
    # 0. Extract training, test data and responses
    TrainSet <- dataShort [trainIndices [[i]], ]
    if(length(trainIndices) != 1) {
      TestSet <- dataShort [-trainIndices [[i]], ]
    }
    else {
      TestSet <- TrainSet
    }

    # 1. Convert training data to long format
    if(!is.null(idColumn)) {
      TrainLong <- dataLongTimeDep(dataSemiLong = TrainSet, timeColumn = as.character(survModelFormula) [2], 
                                    eventColumn = as.character(censModelFormula) [2], idColumn = idColumn,
                                    timeAsFactor = timeAsFactor)
    }
    else {
      TrainLong <- dataLong(dataShort = TrainSet, timeColumn = as.character(survModelFormula) [2], 
                             eventColumn = as.character(censModelFormula) [2], timeAsFactor = timeAsFactor)
    }
    
    # 2. Convert response in training data to censoring variable
    TrainLong <- dataCensoring(dataShort = TrainLong, timeColumn = "timeInt", shortFormat = FALSE)
    
    # 3. Convert test data to long format
    if(!is.null(idColumn)) {
      TestLong <- dataLongTimeDep(dataSemiLong = TestSet, timeColumn = as.character(survModelFormula) [2], 
                                   eventColumn = as.character(censModelFormula) [2], idColumn = idColumn,
                                   timeAsFactor = timeAsFactor)
    } else {
      TestLong <- dataLong(dataShort = TestSet, timeColumn = as.character(survModelFormula) [2], 
                            eventColumn = as.character(censModelFormula) [2], timeAsFactor = timeAsFactor)
    }

    # 4. Fit censoring model on training data in long format
    CensnewFormula <- update(censModelFormula, yCens ~ timeInt + .)
    CensFit <- glm (formula = CensnewFormula, data = TrainLong, family = binomial(link = linkFunc), control = glm.control(maxit = 2500))
    
    # 5. Estimate censoring curves on test data
    # Exclude cases with new factor levels in test data in long format
    Check <- "error" %in% class(tryCatch(predict(CensFit, TestLong, type="response"), error = function (e) e))
    if(Check) {
      
      # Which columns are factors in test data?
      IndexFactor <- which(sapply(1:dim(TestLong)[2], function (x) is.factor(TestLong [, x]))==TRUE)
      
      # What are the levels of these factors?
      TestLevelsFactor <- sapply(IndexFactor, function (x) levels(TestLong [, x]))
      
      # First column does not count (censoring process)
      # What are the levels of the corresponding factors in the training data?
      TrainLevelsFactor <- sapply(IndexFactor, function (x) levels(TrainLong [, x + 1]))
      
      # Which levels of the test data exists in the training data factors?
      InLevelsFactor <- lapply(1:length(TestLevelsFactor), function (x) which((TestLevelsFactor [[x]] %in% TrainLevelsFactor [[x]]) == FALSE))
      ExcludeRows <- lapply (1:length(IndexFactor), function (j) which(TestLong [, IndexFactor [j]] %in% TestLevelsFactor [[j]] [InLevelsFactor [[j]]]))
      ExcludeRows <- do.call(c, ExcludeRows)
      
      # Convert Indices of left out test data to complete data set in long format
      ExcludeRowsConv <- vector("integer", length(ExcludeRows))
      for(j in 1:length(ExcludeRows)) {
        I1 <- sapply(1:dim(TrainLongFull) [1], function(x) TrainLongFull [x, -1] == TestLong [ExcludeRows [j], -1])
        ExcludeRowsConv [j] <- which(sapply(1:dim(I1) [2], function (x) all(I1 [, x])) == TRUE)
      }
      ExcludeRowsCensList [[i]] <- ExcludeRowsConv
      
      # Exclude rows of TestLong
      TestLong <- TestLong [-ExcludeRows, ]
    }
    oneMinuslambdaList [[i]] <- 1 - predict(CensFit, TestLong, type = "response")
    
    # 7. Estimate marker values
    SurvnewFormula <- update(survModelFormula, y ~ timeInt + .)
    SurvFit <- glm (formula = SurvnewFormula, data = TrainLong, family = binomial(link = linkFunc), control = glm.control(maxit = 2500))
    if(timeAsFactor) {
      TestSetExt <- cbind(TestSet, timeInt = factor(TestSet [, as.character(survModelFormula) [2] ]))
      TrainSetExt <- cbind(TrainSet, timeInt = factor(TrainSet [, as.character(survModelFormula) [2] ]))
    }
    else{
      TestSetExt <- cbind(TestSet, timeInt = TestSet [, as.character(survModelFormula) [2] ])
      TrainSetExt <- cbind(TrainSet, timeInt = TrainSet [, as.character(survModelFormula) [2] ])
    }
    
    # Exclude cases with new factor levels in test data in short format
    Check <- "error" %in% class(tryCatch(predict(SurvFit, TestSetExt), error = function (e) e))
    if(Check) {
      # Which columns are factors in test data?
      IndexFactor <- which(sapply(1:dim(TestSetExt)[2], function (x) is.factor(TestSetExt [, x])) == TRUE)
      
      # What are the levels of these factors?
      TestLevelsFactor <- sapply(IndexFactor, function (x) levels(TestSetExt [, x]))
      
      # First column does not count (censoring process)
      # What are the levels of the corresponding factors in the training data?
      TrainLevelsFactor <- sapply(IndexFactor, function (x) levels(TrainSetExt [, x]))
      
      # Which levels of the test data exists in the training data factors?
      InLevelsFactor <- lapply(1:length(TestLevelsFactor), function (x) which((TestLevelsFactor [[x]] %in% TrainLevelsFactor [[x]]) == FALSE))
      ExcludeRows <- lapply (1:length(IndexFactor), function (j) which(TestSetExt [, IndexFactor [j]] %in% TestLevelsFactor [[j]] [InLevelsFactor [[j]] ]))
      ExcludeRows <- do.call(c, ExcludeRows)
      
      # Convert excluded rows of test data in short format to complete data in short format (necessary for newEvent and newTime)
      ExcludeRowsConvShort <- vector("integer", length(ExcludeRows))
      for(j in 1:length(ExcludeRows)) {
        I1 <- sapply(1:dim(dataShort) [1], function(x) dataShort [x, ] == TestSetExt [ExcludeRows [j], -dim(TestSetExt) [2] ])
        ExcludeRowsConvShort [j] <- which(sapply(1:dim(I1) [2], function (x) all(I1 [, x])) == TRUE)
      }
      ExcludeRowsdataShortList [[i]] <- ExcludeRowsConvShort
      
      # Exclude rows of test data in short format
      TestSetExt <- TestSetExt [-ExcludeRows, ]
    }
    markerList [[i]] <- predict(SurvFit, TestSetExt)
  }

  # 8. Estimate sensitivity
  # Estimate survival function of censoring process (complete data set)
  oneMinuslambda <- do.call(c, oneMinuslambdaList)
  # Merge excluded rows
  ExcludeRowsCens <- do.call(c, ExcludeRowsCensList)
  ExcludeRowsdataShort <- do.call(c, ExcludeRowsdataShortList)
  if(!is.null(ExcludeRowsCens)) {
    TrainLongFullExc <- TrainLongFull [-ExcludeRowsCens, ]
  }
  else {
    TrainLongFullExc <- TrainLongFull
  }
  G <- aggregate(oneMinuslambda ~ obj, FUN = cumprod, data = TrainLongFullExc, simplify = FALSE)
  if(!is.null(ExcludeRowsdataShort)) {
    newEvent <- dataShort [-ExcludeRowsdataShort, as.character(censModelFormula) [2]]
    newTime <- dataShort [-ExcludeRowsdataShort, as.character(survModelFormula) [2]]
  }
  else {
    newEvent <- dataShort [, as.character(censModelFormula) [2]]
    newTime <- dataShort [, as.character(survModelFormula) [2]]
  }
  n <- length(newEvent)
  if(is.null(idColumn)) {
    GT <- sapply(1:n, function(u){
      if (newTime[u] > 1)
        return(G[[2]] [u] [[1]] [newTime[u]-1]) else
          return(1) } )
  }
  else{
    GT <- sapply(1:n, function(u){
      if (newTime[u] > 1)
        return(G[[2]] [TrainLongFullExc [dataShort[u, idColumn], "obj"] ] [[1]] [newTime[u] - 1]) else
          return(1) } )
  }
  # Merge markers
  marker <- do.call(c, markerList)
  RET <- sapply(marker, sens)
  orderMarker <- order(marker)
  tempDat <- data.frame(cutoff = marker[orderMarker], tpr = RET[orderMarker])
  rownames(tempDat) <- 1:dim(tempDat) [1]
  
  RET <- list(Output = tempDat, 
              Input = list(timepoint = timepoint, dataShort = dataShort, trainIndices = trainIndices, 
                         survModelFormula = survModelFormula, censModelFormula = censModelFormula, 
                         linkFunc = linkFunc, idColumn = idColumn, Short = FALSE, 
                         timeAsFactor = timeAsFactor, orderMarker = orderMarker))
  class(RET) <- "discSurvTprUno"
  return(RET)
}

#' @rdname tprUno
#' @param x Object of class "discSurvTprUno" ("class discSurvTprUno")
#' @param \dots Additional arguments to the print function
# #' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @keywords internal
#' @method print discSurvTprUno
#' @noRd
print.discSurvTprUno <- function (x, ...) {
  x$Output[, "cutoff"] <- round(x$Output[, "cutoff"], 4)
  if(!any(is.na(x$Output[, "tpr"]))) {
    x$Output[, "tpr"] <- round(x$Output[, "tpr"], 4)
  } 
  print(x$Output, ...)
}

#' @rdname tprUno
#' @param x Object of class "discSurvTprUno" ("class discSurvTprUno")
#' @param \dots Additional arguments to the plot function
# #' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @keywords internal
#' @method plot discSurvTprUno
#' @noRd
plot.discSurvTprUno <- function (x, ...) {
  if(any(is.na(x$Output [, "tpr"]))) {
    return("No plot available, because there are missing values in tpr!")
  }
  plot(x = x$Output [, "cutoff"], y = x$Output [, "tpr"], xlab = "Cutoff", ylab = "Tpr", las = 1, type = "l", main = paste("Tpr(c, t=", x$Input$timepoint, ")", sep = ""), ...)
}

###########################
# tprUnoShort
###########################

# Description
# Computes tprUno given marker values for a general model without implicit estimation

# Input
# timepoint: Timepoint (integer scalar)
# marker: Linear predictor values of each person in the test data
# newTime: Discrete time of each person in the test data
# newEvent: Event indicator of each person in the test data
# trainTime: 
# trainEvent: 

# Output
# data.frame with columns:
  # cutoff:
  # tpr: True positive rate (numeric) \in [0, 1]



#' True Positive Rate for arbitrary predition models
#' 
#' Estimates the true positive rate (based on concept of Uno, et al.) for an
#' arbitrary discrete survival prediction model on one test data set.
#' 
#' This function is useful, if other models than generalized, linear models
#' (glm) should be used for prediction. In the case of glm better use the cross
#' validation version \code{\link{tprUno}}.
#' 
#' @param timepoint Gives the discrete time interval of which the tpr is
#' evaluated ("numeric vector").
#' @param marker Gives the predicted values of the linear predictor of a
#' regression model ("numeric vector"). May also be on the response scale.
#' @param testTime Time intervals in the test data ("integer vector").
#' @param testEvent Event indicators in the test data ("binary vector").
#' @param trainTime Time intervals in the training data ("integer vector").
#' @param trainEvent Event indicators in the training data ("binary vector").
#' @return List with objects \itemize{ \item{Output} Data frame with two columns:
#' \emph{cutoff} gives the different marker values and \emph{fpr} the false positive
#' rates \item{Input} A list of given argument input values (saved for
#' reference). Another list element is \code{selectInd}, which gives the
#' selected indices of the marker values with time intervals available in both
#' training and test sets. In addition there is the list element
#' \code{orderMarker}, which gives the indices of the marker values in
#' increasing order. }
#' @note It is assumed that all time points up to the last observed interval
#' [a_{q-1}, a_q) are available.
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' 
#' Matthias Schmid \email{matthias.schmid@@imbie.uni-bonn.de}
#' @seealso \code{\link[caret]{createDataPartition}}, \code{\link{glm}}
#' @references 
#' \insertRef{schmidDiscMeasure}{discSurv} \cr\cr
#' \insertRef{unoEvalPred}{discSurv} \cr\cr
#' \insertRef{heagertySurvROC}{discSurv}
#' @keywords internal
#' @examples
#' 
#' ##################################################
#' # Example with unemployment data and prior fitting
#' 
#' library(Ecdat)
#' library(caret)
#' library(mgcv)
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
#' UnempDurSubsetTrainLong <- dataLong(dataShort=UnempDurSubsetTrain, 
#' timeColumn="spell", eventColumn="censor1")
#' 
#' # Estimate gam with smooth baseline
#' gamFit <- gam(formula = y ~ s(I(as.numeric(as.character(timeInt)))) + 
#' s(age) + s(logwage), data = UnempDurSubsetTrainLong, family = binomial())
#' gamFitPreds <- predict(gamFit, newdata = cbind(UnempDurSubsetTest, 
#' timeInt = UnempDurSubsetTest$spell))
#' 
#' # Estimate tpr given one training and one test sample
#' tprGamFit <- tprUnoShort (timepoint = 1, marker = gamFitPreds, 
#' testTime = UnempDurSubsetTest$spell, testEvent = UnempDurSubsetTest$censor1, 
#' trainTime = UnempDurSubsetTrain$spell, trainEvent = UnempDurSubsetTrain$censor1)
#' plot(tprGamFit)
#' 
#' #####################################
#' # Example National Wilm's Tumor Study
#' 
#' library(survival)
#' head(nwtco)
#' summary(nwtco$rel)
#' 
#' # Select subset
#' set.seed(-375)
#' Indices <- sample(1:dim(nwtco)[1], 500)
#' nwtcoSub <- nwtco [Indices, ]
#' 
#' # Convert time range to 30 intervals
#' intLim <- quantile(nwtcoSub$edrel, prob = seq(0, 1, length.out = 30))
#' intLim [length(intLim)] <- intLim [length(intLim)] + 1
#' nwtcoSubTemp <- contToDisc(dataShort = nwtcoSub, timeColumn = "edrel", intervalLimits = intLim)
#' nwtcoSubTemp$instit <- factor(nwtcoSubTemp$instit)
#' nwtcoSubTemp$histol <- factor(nwtcoSubTemp$histol)
#' nwtcoSubTemp$stage <- factor(nwtcoSubTemp$stage)
#' 
#' # Split in training and test sample
#' set.seed(-570)
#' TrainingSample <- sample(1:dim(nwtcoSubTemp)[1], round(dim(nwtcoSubTemp)[1]*0.75))
#' nwtcoSubTempTrain <- nwtcoSubTemp [TrainingSample, ]
#' nwtcoSubTempTest <- nwtcoSubTemp [-TrainingSample, ]
#' 
#' # Convert to long format
#' nwtcoSubTempTrainLong <- dataLong(dataShort = nwtcoSubTempTrain, 
#' timeColumn = "timeDisc", eventColumn = "rel")
#' 
#' # Estimate glm
#' inputFormula <- y ~ timeInt + histol + instit + stage
#' glmFit <- glm(formula = inputFormula, data = nwtcoSubTempTrainLong, family = binomial())
#' linPreds <- predict(glmFit, newdata = cbind(nwtcoSubTempTest, 
#' timeInt = nwtcoSubTempTest$timeDisc))
#' 
#' # Estimate tpr given one training and one test sample at time interval 5
#' tprFit <- tprUnoShort(timepoint = 5, marker = linPreds, 
#' testTime = nwtcoSubTempTest$timeDisc, testEvent = nwtcoSubTempTest$rel, 
#' trainTime = nwtcoSubTempTrain$timeDisc, trainEvent = nwtcoSubTempTrain$rel)
#' plot(tprFit)
#' 
#' @noRd
tprUnoShort <- function (timepoint, marker, testTime, testEvent, trainTime, trainEvent) {

  # Construct short data format
  dataSetShort <- data.frame(trainTime = trainTime, trainEvent = trainEvent)
  
  # # Expand training data in long format with censoring variable
  # dataSetLong <- dataLong(dataShort = dataSetShort, timeColumn = "trainTime", 
  #                          eventColumn = "trainEvent")
  # dataSetLongCens <- dataCensoring(dataShort = dataSetLong, respColumn = "y", 
  #                                   timeColumn = "timeInt")

  # Convert back to short format
  dataSetLongCensTrans <- dataCensoring(dataShort = dataSetShort, 
                                              eventColumns = "trainEvent", 
                                              timeColumn = "trainTime")

  # Exclude test time intervals not observed in training data
  markerInput <- marker
  testEventInput <- testEvent
  testTimeInput <- testTime
  selectInd <- testTime %in% intersect(testTime, trainTime)
  marker <- marker[selectInd]
  testEvent <- testEvent[selectInd]
  testTime <- testTime[selectInd]
  
  if(length(testTime) == 0){
    orderMarker <- order(marker)
    tempDat <- data.frame(cutoff = marker[orderMarker], tpr = NA)
    rownames(tempDat) <- 1:dim(tempDat) [1]
    RET <- list(Output = tempDat, Input = list(timepoint = timepoint, marker = markerInput,
                                           testTime = testTimeInput, testEvent = testEventInput,
                                           trainTime = trainTime, trainEvent = trainEvent,
                                           Short = TRUE, selectInd = selectInd,
                                           orderMarker = orderMarker))
    class(RET) <- "discSurvTprUno"
    return(RET)
  }

  # Estimate nonparametric survival function of censoring variable
  tempLifeTab <- lifeTable(dataShort = dataSetLongCensTrans, timeColumn = "timeCens", 
                            eventColumn = "yCens")
  preG <- tempLifeTab [[1]] [, "S"]
  GT <- c(1, preG)
  GT <- GT [testTime]

  # Help function
  sens <- function(k) {
    sensNum <- sum( (marker > k) * (testTime == timepoint) * testEvent / GT)
    sensDenom <- sum( (testTime == timepoint) * testEvent / GT)

    if (sensDenom > 0) {
      return(sensNum / sensDenom)
    } else{
      return(NA)
    }
  }

  RET <- sapply(marker, sens)
  orderMarker <- order(marker)
  tempDat <- data.frame(cutoff = marker[orderMarker], tpr = RET[orderMarker])
  rownames(tempDat) <- 1:dim(tempDat) [1]
  RET <- list(Output = tempDat, Input = list(timepoint = timepoint, marker = markerInput,
                                         testTime = testTimeInput, testEvent = testEventInput,
                                         trainTime = trainTime, trainEvent = trainEvent,
                                         Short = TRUE, selectInd = selectInd,
                                         orderMarker = orderMarker))
  class(RET) <- "discSurvTprUno"
  return(RET)
}

# Old, slow version
# tprUnoShort <- function (timepoint, marker, newTime, newEvent, trainTime, trainEvent) {
#   
#   # Help function
#   TransformLongToShort <- function (dataSetLong, idColumn) {
#     SplitLongData <- split(dataSetLong, dataSetLong [, idColumn])
#     NewDataSplit <- list()
#     for(i in 1:length(SplitLongData)) {
#       if(is.na(tail(SplitLongData [[i]], n = 1) [,"yCens"])) {
#         NewDataSplit [[i]] <- tail(SplitLongData [[i]], n = 2) [-2,]
#       }
#       else {
#         NewDataSplit [[i]] <- tail(SplitLongData [[i]], n = 1)
#       }
#     }
#     result <- do.call(rbind, NewDataSplit)
#     return(result)
#   }
#   
#   # Expand training data in long format with censoring variable
#   dataSetLong <- dataLong(dataShort = data.frame(trainTime = trainTime, trainEvent = trainEvent), timeColumn = "trainTime", eventColumn = "trainEvent")
#   dataSetLongCens <- dataCensoring(dataShort = dataSetLong, respColumn = "y", idColumn = "obj")
# 
#   # Convert back to short format
#   dataSetLongCensTrans <- TransformLongToShort (dataSetLong = dataSetLongCens, idColumn = "obj")
#   dataSetLongCensTrans <- na.omit(dataSetLongCensTrans)
#   
#   # Exclude test time intervals not observed in training data
#   newTimeInput <- newTime
#   newTime <- newTime[newTime %in% intersect(newTime, trainTime)]
#   if(length(newTime) == 0){
#     tempDat <- data.frame(cutoff = sort(marker), tpr = NA)
#     rownames(tempDat) <- 1:dim(tempDat) [1]
#     RET <- list(Output = tempDat, Input = list(timepoint = timepoint, marker = marker, 
#                                                newTime = newTimeInput, newEvent = newEvent, 
#                                                trainTime = trainTime, trainEvent = trainEvent, 
#                                                Short = TRUE))
#     class(RET) <- "discSurvTprUno"
#     return(RET)
#   }
#   
#   # Estimate nonparametric survival function of censoring variable 
#   tempLifeTab <- lifeTable(dataShort = dataSetLongCensTrans, timeColumn = "timeInt", eventColumn = "yCens")
#   preG <- tempLifeTab [[1]] [, "S"]
#   GT <- c(1, preG)
#   GT <- GT [newTime]
#   
#   # Help function
#   sens <- function(k) {
#     sensNum <- sum( (marker > k) * (newTime == timepoint) * newEvent / GT)
#     sensDenom <- sum( (newTime == timepoint) * newEvent / GT)
#     
#     if (sensDenom > 0) {
#       return(sensNum / sensDenom)
#     } else{
#       return(NA)
#     }
#   }
#   
#   RET <- sapply(marker, sens)
#   tempDat <- data.frame(cutoff = sort(marker), tpr = RET [order(marker)])
#   rownames(tempDat) <- 1:dim(tempDat) [1]
#   RET <- list(Output = tempDat, Input = list(timepoint = timepoint, marker = marker, 
#                                          newTime = newTimeInput, newEvent = newEvent, 
#                                          trainTime = trainTime, trainEvent = trainEvent, 
#                                          Short = TRUE))
#   class(RET) <- "discSurvTprUno"
#   return(RET)
# }

#############################
# fprUno
#############################

##############
# Description
# Estimates the predictive false positive rate (fpr) based on cross validation and generalized, linear models

#######
# Input
# timepoint: Discrete time interval given that the false positive rate is evaluated (integer scalar)
# dataSet: Original data. Should be in format data.frame()
# trainIndices: List of Indices from original data used for training (list of integer vectors). 
  # The length of the list is equal to the number of cross valdiation samples
# survModelFormula: Formula of the survival model
# censModelFormula: Formula of the censoring model. Normally this is done without covariates
# linkFunc: Link function of the generalized, linear model see glm
# idColumn: Name of the column with identification numbers of persons. 
# Default NULL means, that each row equals one person (no repeated measurements).

# Output
# data.frame with columns
  # cutoff: Cut off values of the linear predictor (numeric vector)
  # fpr: False positive rate (numeric vector)


#' @name fprUno
#' @title False Positive Rate Uno
#' 
#' @description Estimates the predictive false positive rate (fpr) based on cross validation
#' and generalized, linear models (see \code{\link{glm}}). The concept was
#' suggested by Uno et al. (see references)
#' 
#' 
#' @param timepoint Discrete time interval given that the false positive rate
#' is evaluated ("integer scalar").
#' @param dataShort Original data in short format ("class data.frame").
#' @param trainIndices List of Indices from original data used for training
#' ("vector list"). The length of the list is equal to the number of
#' cross validation samples.
#' @param survModelFormula Formula of the discrete survival model("class formula"). It is used
#' in a generalized, linear model.
#' @param censModelFormula Formula of the censoring model("class formula"). It is used in a
#' generalized, linear model. Usually this is done without covariates.
#' @param linkFunc Link function of the generalized, linear model("character vector").
#' @param idColumn Name of the column with identification numbers of persons("character vector").
#' Default NULL means, that each row equals one person (no repeated
#' measurements).
#' @param timeAsFactor Should the time intervals be coded as factor("logical vector")? Default is
#' to use factor. If the argument is false, the column is coded as numeric.
#' @return \itemize{ \item{Output} List with objects: \itemize{ \item{Output}
#' Data frame with two columns: \emph{cutoff} gives the different marker values and
#' \emph{fpr} the false positive rates \item{Input} A list of given argument input
#' values (saved for reference). In addition there is the list element
#' \code{orderMarker}, which gives the indices of the marker values in
#' increasing order.  } }
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' 
#' Matthias Schmid \email{matthias.schmid@@imbie.uni-bonn.de}
#' @seealso \code{\link[caret]{createDataPartition}}, \code{\link{glm}}
#' @references 
#' \insertRef{schmidDiscMeasure}{discSurv} \cr\cr
#' \insertRef{unoEvalPred}{discSurv} \cr\cr
#' \insertRef{heagertySurvROC}{discSurv}
#' @keywords internal
#' @examples
#' 
#' # Example with cross validation and unemployment data 
#' library(Ecdat)
#' library(caret)
#' data(UnempDur)
#' summary(UnempDur$spell)
#' 
#' # Extract subset of data
#' set.seed(635)
#' IDsample <- sample(1:dim(UnempDur)[1], 100)
#' UnempDurSubset <- UnempDur [IDsample, ]
#' head(UnempDurSubset)
#' range(UnempDurSubset$spell)
#' set.seed(7550)
#' CVfolds <- createFolds (y = UnempDurSubset$spell, returnTrain = TRUE, k = 2)
#' 
#' # Estimate false positive rate of time interval 7:
#' tryFPR <- fprUno (timepoint = 7, dataShort = UnempDurSubset, trainIndices = CVfolds,  
#' survModelFormula = spell ~ age + logwage, censModelFormula = censor1 ~ 1, 
#' linkFunc = "logit", idColumn = NULL, timeAsFactor = FALSE)
#' tryFPR
#' plot(tryFPR)
#' 
#' @noRd
fprUno <- function(timepoint, dataShort, trainIndices, survModelFormula, censModelFormula, linkFunc = "logit", idColumn = NULL, timeAsFactor = TRUE) {
  
  # Input Checks
  if(length(timepoint) != 1 || !(timepoint == floor(timepoint))) {stop("Argument *timepoint* is not in the correct format! Please specify as integer scalar value.")}
  if(!is.data.frame(dataShort)) {stop("Argument *dataShort* is not in the correct format! Please specify as data.frame object.")}
  if(!is.list(trainIndices)) {stop("Argument *trainIndices* is not in the correct format! Please specify a list.")}
  InputCheck1 <- all(sapply(1:length(trainIndices), function (x) is.integer(trainIndices [[x]])))
  if(!InputCheck1) {stop("Sublists of *trainIndices* are not all integer values! Please specify a list of integer Indices.")}
  # Checks if union of test data indices equals all indices in the complete data
  if(length(trainIndices)!=1) {
    InputCheck2 <- all(sort(as.numeric(do.call(c, lapply(trainIndices, function (x) setdiff(1:dim(dataShort) [1], x))))) == (1:dim(dataShort) [1]))
  }
  else {
    InputCheck2 <- all(trainIndices [[1]] == (1:dim(dataShort) [1]))
  }
  if(!InputCheck2) {stop("Argument *trainIndices* does not contain cross validation samples! Please ensure that the union of all test indices equals the indices of the complete data set.")}
  # Formula checks
  if(!("formula" %in% class(censModelFormula))) {stop("*censModelFormula* is not of class formula! Please specify a valid formula, e. g. yCens ~ 1.")}
  if(!("formula" %in% class(survModelFormula))) {stop("*survModelFormula* is not of class formula! Please specify a valid formula, e. g. y ~ x")}
  if(!(any(names(dataShort) == idColumn) | is.null(idColumn))) {stop("Argument *idColumn* is not available in *dataShort*! Please specify the correct column name of the identification number.")}
  
  # Help function
  spec <- function(k){

    specNum <- sum( (marker <= k) * (newTime > timepoint), na.rm = TRUE)
    specDenom <- sum(newTime > timepoint, na.rm = TRUE)

    if (specDenom > 0)
      return(specNum / specDenom) else
        return(0)
  }

  # Loop across all training data sets
  RET <- vector("list", length(trainIndices))
  markerList <- vector("list", length(trainIndices))
  ExcludeRowsdataShortList <- vector("list", length(trainIndices))
  for(i in 1:length(trainIndices)) {
    
    # 0. Extract training, test data and responses
    TrainSet <- dataShort [trainIndices [[i]], ]
    if(length(trainIndices) != 1) {
      TestSet <- dataShort [-trainIndices [[i]], ]
    }
    else {
      TestSet <- TrainSet
    }
    
    # 1. Convert training data to long format
    if(!is.null(idColumn)) {
      TrainLong <- dataLongTimeDep(dataSemiLong = TrainSet, timeColumn = as.character(survModelFormula) [2], 
                                    eventColumn = as.character(censModelFormula) [2], idColumn = idColumn, 
                                    timeAsFactor = timeAsFactor)
    }
    else {
      TrainLong <- dataLong(dataShort = TrainSet, timeColumn = as.character(survModelFormula) [2], 
                             eventColumn = as.character(censModelFormula) [2], timeAsFactor = timeAsFactor)
    }
    
    # 2. Convert response in training data to censoring variable
    TrainLong <- dataCensoring(dataShort = TrainLong, timeColumn = "timeInt", shortFormat = FALSE)
    
    # 3. Convert test data to long format
    if(!is.null(idColumn)) {
      TestLong <- dataLongTimeDep(dataSemiLong = TestSet, timeColumn = as.character(survModelFormula) [2], 
                                   eventColumn = as.character(censModelFormula) [2], idColumn = idColumn, timeAsFactor = timeAsFactor)
    }
    else {
      TestLong <- dataLong(dataShort = TestSet, timeColumn = as.character(survModelFormula) [2], 
                            eventColumn = as.character(censModelFormula) [2], timeAsFactor = timeAsFactor)
    }

    # 7. Estimate marker values
    SurvnewFormula <- update(survModelFormula, y ~ timeInt + .)
    SurvFit <- glm (formula = SurvnewFormula, data = TrainLong, family = binomial(link = linkFunc), control = glm.control(maxit = 2500))
    if(timeAsFactor) {
      TestSetExt <- cbind(TestSet, timeInt = factor(TestSet [, as.character(survModelFormula) [2] ]))
      TrainSetExt <- cbind(TrainSet, timeInt = factor(TrainSet [, as.character(survModelFormula) [2] ]))
    }
    else{
      TestSetExt <- cbind(TestSet, timeInt = TestSet [, as.character(survModelFormula) [2] ])
      TrainSetExt <- cbind(TrainSet, timeInt = TrainSet [, as.character(survModelFormula) [2] ])
    }

    Check <- "error" %in% class(tryCatch(predict(SurvFit, TestSetExt), error = function (e) e))
    if(Check) {
      # Which columns are factors in test data?
      IndexFactor <- which(sapply(1:dim(TestSetExt)[2], function (x) is.factor(TestSetExt [, x]))==TRUE)
      
      # What are the levels of these factors?
      TestLevelsFactor <- sapply(IndexFactor, function (x) levels(TestSetExt [, x]))
      
      # First column does not count (censoring process)
      # What are the levels of the corresponding factors in the training data?
      TrainLevelsFactor <- sapply(IndexFactor, function (x) levels(TrainSet [, x]))
      
      # Which levels of the test data exists in the training data factors?
      InLevelsFactor <- lapply(1:length(TestLevelsFactor), function (x) which((TestLevelsFactor [[x]] %in% TrainLevelsFactor [[x]])==FALSE))
      ExcludeRows <- lapply (1:length(IndexFactor), function (j) which(TestSetExt [, IndexFactor [j]] %in% TestLevelsFactor [[j]] [InLevelsFactor [[j]] ]))
      ExcludeRows <- do.call(c, ExcludeRows)
      
      # Convert excluded rows of test data in short format to complete data in short format (necessary for newEvent and newTime)
      ExcludeRowsConvShort <- vector("integer", length(ExcludeRows))
      for(j in 1:length(ExcludeRows)) {
        I1 <- sapply(1:dim(dataShort) [1], function(x) dataShort [x, ] == TestSetExt [ExcludeRows [j], - dim(TestSetExt) [2] ])
        ExcludeRowsConvShort [j] <- which(sapply(1:dim(I1) [2], function (x) all(I1 [, x]))==TRUE)
      }
      ExcludeRowsdataShortList [[i]] <- ExcludeRowsConvShort
      
      TestSetExt <- TestSetExt [-ExcludeRows, ]
    }
    markerList [[i]] <- predict(SurvFit, TestSetExt)
  }
  
  # 8. Estimate specificity
  marker <- do.call(c, markerList)
  ExcludeRowsdataShort <- do.call(c, ExcludeRowsdataShortList)
  if(!is.null(ExcludeRowsdataShort)) {
    newEvent <- dataShort [-ExcludeRowsdataShort, as.character(censModelFormula) [2]]
    newTime <- dataShort [-ExcludeRowsdataShort, as.character(survModelFormula) [2]]
  }
  else {
    newEvent <- dataShort [, as.character(censModelFormula) [2]]
    newTime <- dataShort [, as.character(survModelFormula) [2]]
  }
  RET <- sapply(marker, spec)
  
  # Output
  orderMarker <- order(marker)
  tempDat <- data.frame(cutoff = marker[orderMarker], fpr = 1 - RET[orderMarker])
  rownames(tempDat) <- 1:dim(tempDat) [1]
  RET <- list(Output = tempDat, 
              Input = list(timepoint = timepoint, dataShort = dataShort, trainIndices = trainIndices, 
                         survModelFormula = survModelFormula, censModelFormula = censModelFormula, 
                         linkFunc = linkFunc, idColumn = idColumn, Short = FALSE, 
                         timeAsFactor = timeAsFactor, orderMarker = orderMarker))
  class(RET) <- "discSurvFprUno"
  return(RET)
}

#' @rdname fprUno
#' @param x Object of class "discSurvFprUno"("class discSurvFprUno")
# #' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @keywords internal
#' @method print discSurvFprUno
#' @noRd
print.discSurvFprUno <- function (x, ...) {
  x$Output[, "cutoff"] <- round(x$Output[, "cutoff"], 4)
  if(!any(is.na(x$Output[, "fpr"]))) {
    x$Output[, "fpr"] <- round(x$Output[, "fpr"], 4)
  } 
  print(x$Output, ...)
}

#' @rdname fprUno
#' @param x Object of class "discSurvFprUno"("class discSurvFprUno")
#' @param \dots Additional arguments to S3 methods.
# #' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @keywords internal
#' @method plot discSurvFprUno
#' @noRd
plot.discSurvFprUno <- function (x, ...) {
  if(any(is.na(x$Output [, "fpr"]))) {
    return("No plot available, because there are missing values in tpr!")
  }
  plot(x = x$Output [, "cutoff"], y = x$Output [, "fpr"], xlab = "Cutoff", ylab = "Fpr", las = 1, type = "l", main = paste("Fpr(c, t=", x$Input$timepoint, ")", sep = ""), ...)
}

#######################
# fprUnoShort
#######################

# Description
# Estimates the false positive rate given prior estimated marker values



#' False Positive Rate for arbitrary predition models
#' 
#' Estimates the false positive rate (based on concept of Uno, et al.) for an
#' arbitrary discrete survival prediction model on one test data set.
#' 
#' This function is useful, if other models than generalized, linear models
#' (glm) should be used for prediction. In the case of glm better use the cross
#' validation version \code{\link{fprUno}}.
#' 
#' @param timepoint Gives the discrete time interval of which the fpr is
#' evaluated ("numeric vector").
#' @param marker Gives the predicted values of the linear predictor of a
#' regression model ("numeric vector"). May also be on the response scale.
#' @param testTime New time intervals in the test data ("integer vector").
#' @return \itemize{ \item{Output} A list with objects: \itemize{ \item{Output}
#' Data frame with two columns: "cutoff" gives the different marker values and
#' \emph{fpr} the false positive rates \item{Input} A list of given argument input
#' values (saved for reference). In addition there is the list element
#' \code{orderMarker}, which gives the indices of the marker values in
#' increasing order.  } }
#' @note It is assumed that all time points up to the last observed interval
#' [a_{q-1}, a_q) are available.
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' 
#' Matthias Schmid \email{matthias.schmid@@imbie.uni-bonn.de}
#' @seealso \code{\link[caret]{createDataPartition}}, \code{\link{glm}}
#' @references 
#' \insertRef{schmidDiscMeasure}{discSurv} \cr\cr
#' \insertRef{unoEvalPred}{discSurv} \cr\cr
#' \insertRef{heagertySurvROC}{discSurv}
#' @keywords internal
#' @examples
#' 
#' ##################################################
#' # Example with unemployment data and prior fitting
#' 
#' library(Ecdat)
#' library(caret)
#' library(mgcv)
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
#' UnempDurSubsetTrainLong <- dataLong(dataShort = UnempDurSubsetTrain, 
#' timeColumn = "spell", eventColumn = "censor1")
#' 
#' # Estimate gam with smooth baseline
#' gamFit <- gam(formula = y ~ s(I(as.numeric(as.character(timeInt)))) + 
#' s(age) + s(logwage), data = UnempDurSubsetTrainLong, family = binomial())
#' gamFitPreds <- predict(gamFit, newdata = cbind(UnempDurSubsetTest, timeInt = UnempDurSubsetTest$spell))
#' 
#' # Estimate false positive rate
#' fprGamFit <- fprUnoShort (timepoint = 1, marker = gamFitPreds, 
#' testTime = UnempDurSubsetTest$spell)
#' plot(fprGamFit)
#' 
#' #####################################
#' # Example National Wilm's Tumor Study
#' 
#' library(survival)
#' head(nwtco)
#' summary(nwtco$rel)
#' 
#' # Select subset
#' set.seed(-375)
#' Indices <- sample(1:dim(nwtco)[1], 500)
#' nwtcoSub <- nwtco [Indices, ]
#' 
#' # Convert time range to 30 intervals
#' intLim <- quantile(nwtcoSub$edrel, prob = seq(0, 1, length.out = 30))
#' intLim [length(intLim)] <- intLim [length(intLim)] + 1
#' nwtcoSubTemp <- contToDisc(dataShort = nwtcoSub, timeColumn = "edrel", intervalLimits = intLim)
#' nwtcoSubTemp$instit <- factor(nwtcoSubTemp$instit)
#' nwtcoSubTemp$histol <- factor(nwtcoSubTemp$histol)
#' nwtcoSubTemp$stage <- factor(nwtcoSubTemp$stage)
#' 
#' # Split in training and test sample
#' set.seed(-570)
#' TrainingSample <- sample(1:dim(nwtcoSubTemp)[1], round(dim(nwtcoSubTemp)[1]*0.75))
#' nwtcoSubTempTrain <- nwtcoSubTemp [TrainingSample, ]
#' nwtcoSubTempTest <- nwtcoSubTemp [-TrainingSample, ]
#' 
#' # Convert to long format
#' nwtcoSubTempTrainLong <- dataLong(dataShort = nwtcoSubTempTrain, 
#' timeColumn = "timeDisc", eventColumn = "rel")
#' 
#' # Estimate glm
#' inputFormula <- y ~ timeInt + histol + instit + stage
#' glmFit <- glm(formula = inputFormula, data = nwtcoSubTempTrainLong, family = binomial())
#' linPreds <- predict(glmFit, newdata = cbind(nwtcoSubTempTest, 
#' timeInt = nwtcoSubTempTest$timeDisc))
#' 
#' # Estimate tpr given one training and one test sample at time interval 10
#' fprFit <- fprUnoShort (timepoint = 10, marker = linPreds, 
#' testTime = nwtcoSubTempTest$timeDisc)
#' plot(fprFit)
#' 
#' @noRd
fprUnoShort <- function (timepoint, marker, testTime) {
  
  # Help function
  spec <- function(k){
    
    specNum <- sum( (marker <= k) * (testTime > timepoint) )
    specDenom <- sum(testTime > timepoint)
    
    if (specDenom > 0) {
      return(specNum / specDenom)
    } else {
      return(NA)
    }
  }
  
  # Output
  RET <- sapply(marker, spec)
  orderMarker <- order(marker)
  tempDat <- data.frame(cutoff = marker[orderMarker], fpr = 1 - RET[orderMarker])
  rownames(tempDat) <- 1:dim(tempDat) [1]
  RET <- list(Output = tempDat, Input = list(timepoint = timepoint, marker = marker, 
                                         testTime = testTime, Short = TRUE,
                                         orderMarker = orderMarker))
  class(RET) <- "discSurvFprUno"
  return(RET)
}

######################
# aucUno
######################

# Description
# Computes the auc (area under the curve) measure given timepoint t
# Both input object must have identical input parameters, e. g. same timepoint, formula!

# Input
# tprObj: Object of class "discSurvTprUno"
# fprObj: Object of class "discSurvFprUno"

# Output
# auc value (numeric scalar) given timepoint of *tprObj* and *fprObj*



#' Area under the Curve Estimation
#' 
#' Estimates the time dependent area under the curve given calculated true
#' positive rate and false positive rate. Both objects ("tprObj", "fprObj") and
#' must have identical input arguments, e. g. same relationship of discrete
#' response and covariates and supplied data sources. The values should be
#' above 0.5 for a well predicting model, because random guessing would get
#' this score.
#' 
#' The auc is estimated by numerical integration of the ROC curve.
#' 
#' @param tprObj Object of class "discSurvTprUno"("class discSurvTprUno"). See function
#' \code{\link{tprUno}}
#' @param fprObj Object of class "discSurvFprUno"("class discSurvFprUno"). See function
#' \code{\link{fprUno}}
#' @return \itemize{ \item{Output} A list with objects: \itemize{ \item{Output}
#' Named numeric vector with auc value \item{Input} A list of given argument
#' input values (saved for reference) } }
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' 
#' Matthias Schmid \email{matthias.schmid@@imbie.uni-bonn.de}
#' @seealso \code{\link{tprUno}}, \code{\link{fprUno}}
#' @references 
#' \insertRef{schmidDiscMeasure}{discSurv} \cr\cr
#' \insertRef{unoEvalPred}{discSurv} \cr\cr
#' \insertRef{heagertySurvROC}{discSurv}
#' @keywords internal
#' @examples
#' 
#' #####################################################
#' # Example with cross validation and unemployment data 
#' 
#' library(Ecdat)
#' library(caret)
#' data(UnempDur)
#' summary(UnempDur$spell)
#' 
#' # Extract subset of data
#' set.seed(635)
#' IDsample <- sample(1:dim(UnempDur)[1], 100)
#' UnempDurSubset <- UnempDur [IDsample, ]
#' head(UnempDurSubset)
#' range(UnempDurSubset$spell)
#' set.seed(7550)
#' CVfolds <- createFolds (y = UnempDurSubset$spell, returnTrain = TRUE, k = 2)
#' 
#' # Estimate true positive rate of time interval 7: 
#' # Correspondes to three and a half month duration (each interval is of length two weeks)
#' tryTPR <- tprUno(timepoint = 7, dataShort = UnempDurSubset, trainIndices = CVfolds, 
#' survModelFormula = spell ~ age + logwage, censModelFormula = censor1 ~ 1, 
#' linkFunc = "logit", idColumn = NULL, timeAsFactor = FALSE)
#' tryTPR
#' plot(tryTPR)
#' 
#' # Estimate false positive rate of time interval 7:
#' tryFPR <- fprUno (timepoint = 7, dataShort = UnempDurSubset, trainIndices = CVfolds,  
#' survModelFormula = spell ~ age + logwage, censModelFormula = censor1 ~ 1, 
#' linkFunc = "logit", idColumn = NULL, timeAsFactor = FALSE)
#' tryFPR
#' plot(tryFPR)
#' 
#' # Estimate auc
#' tryAUC <- aucUno (tprObj = tryTPR, fprObj = tryFPR)
#' tryAUC
#' plot(tryAUC)
#' 
#' #####################################
#' # Example National Wilm's Tumor Study
#' 
#' library(survival)
#' head(nwtco)
#' summary(nwtco$rel)
#' 
#' # Select subset
#' set.seed(-375)
#' Indices <- sample(1:dim(nwtco)[1], 500)
#' nwtcoSub <- nwtco [Indices, ]
#' 
#' # Convert time range to 30 intervals
#' intLim <- quantile(nwtcoSub$edrel, prob = seq(0, 1, length.out = 30))
#' intLim [length(intLim)] <- intLim [length(intLim)] + 1
#' nwtcoSubTemp <- contToDisc(dataShort = nwtcoSub, timeColumn = "edrel", intervalLimits = intLim)
#' nwtcoSubTemp$instit <- factor(nwtcoSubTemp$instit)
#' nwtcoSubTemp$histol <- factor(nwtcoSubTemp$histol)
#' nwtcoSubTemp$stage <- factor(nwtcoSubTemp$stage)
#' 
#' # Split in training and test sample
#' set.seed(-570)
#' TrainingSample <- sample(1:dim(nwtcoSubTemp)[1], round(dim(nwtcoSubTemp)[1]*0.75))
#' nwtcoSubTempTrain <- nwtcoSubTemp [TrainingSample, ]
#' nwtcoSubTempTest <- nwtcoSubTemp [-TrainingSample, ]
#' 
#' # Convert to long format
#' nwtcoSubTempTrainLong <- dataLong(dataShort = nwtcoSubTempTrain, 
#' timeColumn = "timeDisc", eventColumn = "rel")
#' 
#' # Estimate glm
#' inputFormula <- y ~ timeInt + histol + instit + stage
#' glmFit <- glm(formula = inputFormula, data = nwtcoSubTempTrainLong, family = binomial())
#' linPreds <- predict(glmFit, newdata = cbind(nwtcoSubTempTest, 
#' timeInt = nwtcoSubTempTest$timeDisc))
#' 
#' # Estimate tpr given one training and one test sample at time interval 5
#' tprFit <- tprUnoShort(timepoint = 5, marker = linPreds, 
#' newTime = nwtcoSubTempTest$timeDisc, newEvent = nwtcoSubTempTest$rel, 
#' trainTime = nwtcoSubTempTrain$timeDisc, trainEvent = nwtcoSubTempTrain$rel)
#' 
#' # Estimate fpr given one training and one test sample at time interval 5
#' fprFit <- fprUnoShort (timepoint = 5, marker = linPreds, 
#' newTime = nwtcoSubTempTest$timeDisc)
#' 
#' # Estimate auc
#' tryAUC <- aucUno (tprObj = tprFit, fprObj = fprFit)
#' tryAUC
#' plot(tryAUC)
#' 
#' @noRd
aucUno <- function (tprObj, fprObj) {
  
  # Input checks
  if(class(tprObj) != "discSurvTprUno") {stop("This object has not the appropriate class! Please specify an object of class *discSurvTprUno*.")}
  if(class(fprObj) != "discSurvFprUno") {stop("This object has not the appropriate class! Please specify an object of class *discSurvFprUno*.")}
  if(tprObj$Input$Short != fprObj$Input$Short) {stop("Tpr and fpr were computed using different functions! Please ensure that both are estimated either by the cross validated version or the short version.")}
  if(!tprObj$Input$Short) {
    InputCheck <- identical(tprObj$Input, fprObj$Input)
    if(!InputCheck) {stop("Some input parameters of *tprObj* or *fprObj* are not identical! Please check if both objects were estimated using exact identical input values.")}
  }
  else {
    InputCheck1 <- identical(tprObj$Input$timepoint, fprObj$Input$timepoint)
    InputCheck2 <- identical(tprObj$Input$marker, fprObj$Input$marker)
    InputCheck <- all(InputCheck1, InputCheck2)
    if(!InputCheck) {stop("Some input parameters of *tprObj* or *fprObj* are not identical! Please check if both objects were estimated using exact identical input values.")}
  }
  
  # If short format, select specificities according 
  if(tprObj$Input$Short){
    tpr <- c(1, tprObj$Output$tpr)
    fpr <- c(1, fprObj$Output$fpr[ tprObj$Input$selectInd[tprObj$Input$orderMarker] ])
  } else{
    tpr <- c(1, tprObj$Output$tpr)
    fpr <- c(1, fprObj$Output$fpr)
  }
  
  trapz <- function (x, y){ # from package caTools
    idx = 2:length(x)
    return(as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2)
  }
  
  Output <- - trapz(fpr, tpr)
  names(Output) <- paste("AUC(t=", tprObj$Input$timepoint, ")", sep="")
  auc <- list(Output = Output, Input = list(tprObj = tprObj, fprObj = fprObj))
  class(auc) <- "discSurvAucUno"
  return(auc)
}

#' @rdname aucUno
#' @param x Object of class "discSurvAdjDevResid"("class discSurvAdjDevResid")
# #' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @keywords internal
#' @method print discSurvAucUno
#' @noRd
print.discSurvAucUno <- function (x, ...) {
  print(round(x$Output, 4))
}

#' @rdname aucUno
#' @param x Object of class "discSurvAucUno"("class discSurvAdjDevResid")
#' @param \dots Additional arguments to S3 methods.
# #' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @keywords internal
#' @method plot discSurvAucUno
#' @noRd
plot.discSurvAucUno <- function (x, ...) {
  
  # In the short case exclude observations with test times not available
  # in training times
  tprVal <- x$Input$tprObj$Output [, "tpr"]
  if(x$Input$tprObj$Input$Short) {
    fprVal <- x$Input$fprObj$Output [ x$Input$tprObj$Input$selectInd[
      x$Input$tprObj$Input$orderMarker], "fpr"]
  } else{
    fprVal <- x$Input$fprObj$Output [, "fpr"]
  }
  
  if(any(is.na(tprVal)) | any(is.na(fprVal))) {
    return("No plot available, because either tprVal or fprVal contains missing values!")
  }
  plot(x = fprVal, y = tprVal, xlab = "Fpr", ylab = "Tpr", las = 1, type = "l", 
       main = paste("ROC(c, t=", x$Input$tprObj$Input$timepoint, ")", sep = ""), ...)
  lines(x = seq(0, 1, length.out = 500), y = seq(0, 1, length.out = 500), lty = 2)
}

########################################
# concordanceIndex
########################################

# Description
# Calculates the Concordance index (independent measure of time)

# Input
# auc: Double vector with length equal to the number of time intervals
# SurvT: Double vector of estimated survival curve P(T>t) for each time interval

# Output
# Weighted integrated auc over time as measure of accuracy. The format is a double scalar value.



#' Concordance Index
#' 
#' Calculates the concordance index for discrete survival models (independent
#' measure of time). This is the probability that, for a pair of randomly
#' chosen comparable samples, the sample with the higher risk prediction will
#' experience an event before the other sample.
#' 
#' The algorithm extracts all necessary information of the auc object (e. g.
#' marginal probabilities and survival functions).
#' 
#' @param aucObj Object of class "discSurvAucUno"("class discSurvAucUno"). This object is created using
#' the function \code{\link{aucUno}}.
#' @param printTimePoints Should messages be printed for each calculation of a
#' discrete time interval? ("logical vector") Default is FALSE.
#' @return List with objects \itemize{ \item{Output} Concordance index (named
#' numeric vector) \item{Input} List with all input arguments (saved for
#' reference) }
#' @note It is assumed that all time points up to the last observed interval
#' [a_{q-1}, a_q) are available. If not already present, these can be added
#' manually.
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' 
#' Matthias Schmid \email{matthias.schmid@@imbie.uni-bonn.de}
#' @references 
#' \insertRef{schmidDiscMeasure}{discSurv} \cr\cr
#' \insertRef{unoEvalPred}{discSurv} \cr\cr
#' \insertRef{heagertySurvROC}{discSurv}
#' @keywords internal
#' @examples
#' 
#' 
#' # Example with cross validation and unemployment data 
#' library(Ecdat)
#' library(caret)
#' data(UnempDur)
#' 
#' # Evaluation of short term prediction for re-employed at full-time job
#' # Last interval q = 14
#' # -> Set all time points with spell > 13 to time interval 13 and censored
#' lastObsInterval <- 13
#' UnempDurSubset <- UnempDur
#' UnempDurSubset[UnempDurSubset$spell > lastObsInterval, "censor1"] <- 0
#' UnempDurSubset[UnempDurSubset$spell > lastObsInterval, "spell"] <- lastObsInterval
#' head(UnempDurSubset)
#' range(UnempDurSubset$spell)
#' 
#' # Select cross-validation samples
#' set.seed(7550)
#' CVfolds <- createFolds (y = UnempDurSubset$spell, returnTrain = TRUE, k = 2)
#' 
#' # Continuation ratio model formula
#' contModForm <- spell ~ logwage + ui + logwage*ui + age
#' 
#' # Estimate true positive rate of time interval 6: 
#' # Correspondes to three and a half month duration (each interval is of length two weeks)
#' tryTPR <- tprUno (timepoint = 6, dataShort = UnempDurSubset, trainIndices = CVfolds, 
#' survModelFormula = contModForm, censModelFormula = censor1 ~ 1, 
#' linkFunc = "logit", idColumn = NULL, timeAsFactor = FALSE)
#' tryTPR
#' plot(tryTPR)
#' 
#' # Estimate false positive rate of time interval 6:
#' tryFPR <- fprUno (timepoint = 6, dataShort = UnempDurSubset, trainIndices = CVfolds,  
#' survModelFormula = contModForm, censModelFormula = censor1 ~ 1, 
#' linkFunc = "logit", idColumn = NULL, timeAsFactor = FALSE)
#' tryFPR
#' plot(tryFPR)
#' 
#' # Estimate AUC rate of time interval 6:
#' tryAUC <- aucUno (tprObj = tryTPR, fprObj = tryFPR)
#' tryAUC
#' plot(tryAUC)
#' 
#' \dontrun{# Estimate global concordance index:
#' tryConcorIndex <- concorIndex (tryAUC, printTimePoints = TRUE)
#' tryConcorIndex
#' summary(tryConcorIndex)}
#' 
#' @noRd
concorIndex <- function (aucObj, printTimePoints = FALSE) {

  # Input checks
  if(class(aucObj) != "discSurvAucUno") {stop("This object has not the appropriate class! Please specify an object of class *discSurvAucUno*.")}

  if(aucObj$Input$tprObj$Input$Short) {
    
    # Estimate AUC for all t
    marker <- aucObj$Input$tprObj$Input$marker
    newTime <- aucObj$Input$tprObj$Input$testTime
    newEvent <- aucObj$Input$tprObj$Input$testEvent
    trainTime <- aucObj$Input$tprObj$Input$trainTime
    trainEvent <- aucObj$Input$tprObj$Input$trainEvent
    MaxTime <- max(trainTime)-1
    AUCalltime <- vector("numeric", MaxTime)
    for(i in 1:MaxTime) {
      tempTPR <- tprUnoShort (timepoint = i, marker = marker, testTime = newTime, testEvent = newEvent, 
                              trainTime = trainTime, trainEvent = trainEvent)
      tempFPR <- fprUnoShort (timepoint = i, marker = marker, testTime = newTime)
      AUCalltime [i] <- as.numeric(aucUno (tprObj = tempTPR, fprObj = tempFPR)$Output)
      if(printTimePoints) {cat("Progress:", round(i/MaxTime*100, 2), "%;", "Timepoint =", i, "\n")}
    }

    # Estimate nonparametric survival function S(T=t) and marginal probabilities P(T=t)
    tempLifeTab <- lifeTable(dataShort = data.frame(trainTime = trainTime, trainEvent = trainEvent), timeColumn = "trainTime", eventColumn = "trainEvent")
    MargHaz <- tempLifeTab [[1]] [, "hazard"]
    MargSurv <- estSurv(MargHaz)
    MargProb <- estMargProb(MargHaz)
  }
  
  else {

    # Estimate AUC for all t
    MaxTime <- max(aucObj$Input$tprObj$Input$dataSet [, 
                                                      as.character(aucObj$Input$tprObj$Input$survModelFormula) [2] ])-1
    DataSet <- aucObj$Input$tprObj$Input$dataSet
    TrainIndices <- aucObj$Input$tprObj$Input$trainIndices
    SurvModelFormula <- aucObj$Input$tprObj$Input$survModelFormula
    CensModelFormula <- aucObj$Input$tprObj$Input$censModelFormula
    LinkFunc <- aucObj$Input$tprObj$Input$linkFunc
    IdColumn <- aucObj$Input$tprObj$Input$idColumn
    timeAsFactor <- aucObj$Input$tprObj$Input$timeAsFactor
    AUCalltime <- vector("numeric", MaxTime)
    for(i in 1:MaxTime) {
      tempTPR <- tprUno(timepoint = i, dataShort = DataSet, 
                         trainIndices = TrainIndices, 
                         survModelFormula = SurvModelFormula, censModelFormula = CensModelFormula, linkFunc = LinkFunc, idColumn = IdColumn, timeAsFactor = timeAsFactor)
      tempFPR <- fprUno(timepoint = i, dataShort = DataSet, 
                         trainIndices = TrainIndices, 
                         survModelFormula = SurvModelFormula, censModelFormula = CensModelFormula, linkFunc = LinkFunc, idColumn = IdColumn, timeAsFactor = timeAsFactor)
      AUCalltime [i] <- as.numeric(aucUno (tprObj = tempTPR, 
                                           fprObj = tempFPR)$Output)
      if(printTimePoints) {cat("Progress:", round(i/MaxTime*100, 2), "%;", "Timepoint =", i, "\n")}
    }
  
    # Estimate survival function and marginal probabilities without covariates
    if(!is.null(IdColumn)) {
      TrainLongFull <- dataLongTimeDep(dataSemiLong = DataSet, 
                                        timeColumn = as.character(
                                          SurvModelFormula) [2], 
                                        eventColumn = as.character(
                                          CensModelFormula) [2], idColumn = IdColumn, timeAsFactor = timeAsFactor)
    }
    else {
      TrainLongFull <- dataLong(dataShort = DataSet, 
                                 timeColumn = as.character(
                                   SurvModelFormula) [2], 
                                 eventColumn = as.character(
                                   CensModelFormula) [2], 
                                 timeAsFactor = timeAsFactor)
    }
    
    # Estimate marginal survival probability with glm
    MargFormula <- y ~ timeInt
    MargFit <- glm (formula = MargFormula, data = TrainLongFull, 
                    family = binomial(link = LinkFunc), 
                    control = glm.control(maxit = 2500))
    if(timeAsFactor) {
      PredMargData <- data.frame(timeInt = factor(
        min(TrainLongFull [, as.character(SurvModelFormula) [2] ]):
          max(TrainLongFull [, as.character(SurvModelFormula) [2] ])))
    }
    else{
      PredMargData <- data.frame(timeInt = min(
        TrainLongFull [, as.character(SurvModelFormula) [2] ]):
          max(TrainLongFull [, as.character(SurvModelFormula) [2] ]))
    }
    MargHaz <- as.numeric(predict(MargFit, PredMargData, type="response"))
    # Survival function S(T = t) = P(T>t) = \prod_{j=1}^{t1} (1-\lambda (T=j))
    MargSurv <- estSurv(MargHaz)
    # Marginal Probability P(T=t) = \lambda (T=t) \prod_{j=1}^{t-1} (1-\lambda (T=j))
    MargProb <- estMargProb(MargHaz)
  }
  
  # Calculate concordance index
  weights1 <- MargProb * MargSurv / sum(MargProb * MargSurv)
  # Last weight is zero and can therefore be omitted
  # Check for NA, NaN, Inf, -Inf in AUCalltime and weights1
  # Exclude values, if one of them are missing
  AUCind <- is.finite(AUCalltime) & is.finite(weights1[-length(weights1)])
  Concor <- sum( AUCalltime[AUCind] * weights1[-length(weights1)][AUCind] )/
    sum( weights1[-length(weights1)][AUCind] )
  names(Concor) <- "C*"
  names(AUCalltime) <- paste("AUC(t=", 1:MaxTime, "|x)", sep = "")
  Output <- list(Output = Concor, Input = list(aucObj = aucObj, AUC = AUCalltime, 
                                           MargProb = MargProb, MargSurv = MargSurv))
  class(Output) <- "discSurvConcorIndex"
  return(Output)
}

#' @rdname concorIndex
#' @param x Object of class "discSurvConcorIndex"("class discSurvConcorIndex")
# #' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @keywords internal
#' @method print discSurvConcorIndex
#' @noRd
print.discSurvConcorIndex <- function (x, ...) {
  print(round(x$Output, 4))
}

#' @rdname concorIndex
#' @param object Object of class "discSurvConcorIndex"("class discSurvConcorIndex")
#' @param \dots Additional arguments to S3 methods
# #' @author Thomas Welchowski
#' @keywords internal
#' @method summary discSurvConcorIndex
#' @noRd
summary.discSurvConcorIndex <- function (object, ...) {
  cat("Concordance: Should be higher than 0.5 (random assignment)", "\n")
  print(round(object$Output, 4))
  cat("\n", "AUC(t): Should be higher than 0.5 for all time points (random assignment)", "\n")
  print(round(object$Input$AUC, 4))
  cat("\n", "Marginal P(T=t) without covariates (used in weighting)", "\n")
  print(round(object$Input$MargProb, 4))
  cat("\n", "Marginal S(T=t) without covariates (used in weighting)", "\n")
  print(round(object$Input$MargSurv, 4))
}

######################
# predErrDisc

# # 1. Estimate survival function of model
# # 2. Estimate censoring survival function of a covariate free model
# # 3. Calculate observed survival function
# # 4. Calculate prediction error curve given timepoint t
# 
# Problems with this function:
# - Survival function for censoring process is only evaluated with training data!
# - 
# - Function is too complex! Should be simplified! 
# 
# predErrDisc <- function(timepoints, dataSet, trainIndices, survModelFormula, censModelFormula, linkFunc="logit", idColumn=NULL) {
#   
#   # Input Checks
#   if(!is.vector(timepoints) | !all(timepoints==floor(timepoints))) {stop("Argument *timepoints* is not in the correct format! Please specify as integer vector.")}
#   if(!is.data.frame(dataSet)) {stop("Argument *dataSet* is not in the correct format! Please specify as data.frame object.")}
#   if(!is.list(trainIndices)) {stop("Argument *trainIndices* is not in the correct format! Please specify a list.")}
#   InputCheck1 <- all(sapply(1:length(trainIndices), function (x) is.integer(trainIndices [[x]])))
#   if(!InputCheck1) {stop("Sublists of *trainIndices* are not all integer values! Please specify a list of integer Indices.")}
#   if(length(trainIndices)!=1) {
#     InputCheck2 <- all(sort(as.numeric(do.call(c, lapply(trainIndices, function (x) setdiff(1:dim(dataSet) [1], x)))))==(1:dim(dataSet) [1]))
#   }
#   else {
#     InputCheck2 <- all(trainIndices [[1]]==(1:dim(dataSet) [1]))
#   }
#   if(!InputCheck2) {stop("Argument *trainIndices* does not contain cross validation samples! Please ensure that the union of all test indices equals the indices of the complete data set.")}
#   if(!("formula" %in% class(censModelFormula))) {stop("*censModelFormula* is not of class formula! Please specify a valid formula, e. g. yCens ~ 1")}
#   if(!("formula" %in% class(survModelFormula))) {stop("*survModelFormula* is not of class formula! Please specify a valid formula, e. g. y ~ x")}
#   if(!(any(names(dataSet)==idColumn) | is.null(idColumn))) {stop("Argument *idColumn* is not available in *dataSet*! Please specify the correct column name of the identification number.")}
#   
#   Help functions
#     ObsSurvFunc <- function (x) {
#       resp <- TrainLongFullExcSplit [[x]]
#       if(resp [length(resp)]==1) {
#         temp <- c(rep(1, length(resp)-1), 0)
#         result <- tryCatch(temp [timepoints [k] ], error=function (e) NA)
#         return(result)
#       }
#       else {
#         temp <- rep(1, length(resp))
#         result <- tryCatch(temp [timepoints [k] ], error=function (e) NA)
#         return (result)
#       }
#     }
#   
#   WeightFunction <- function () {
#     PartialSum1 <- newEvent * (1 - Sobs) / GT
#     PartialSum2 <- Sobs / GTfixed
#     return(PartialSum1 + PartialSum2)
#   }
#   
#   predErr <- function () {
#     sum(weights [[k]] * (STfixed - Sobs)^2, na.rm=TRUE) / length(weights [[k]] [!is.na(weights [[k]])])
#   }
#   
#   # Loop across all training data sets
#   ExcludeRowsCensList <- vector("list", length(trainIndices))
#   ExcludeRowsDataSetList <- vector("list", length(trainIndices))
#   oneMinusLambdaList <- vector("list", length(trainIndices))
#   oneMinusLambdaSurvList <- vector("list", length(trainIndices))
#   
#   # Convert full sample to long format
#   if(!is.null(idColumn)) {
#     TrainLongFull <- dataLongTimeDep(dataSemiLong=dataSet, timeColumn=as.character(survModelFormula) [2], eventColumn=as.character(censModelFormula) [2], idColumn=idColumn)
#   }
#   else {
#     TrainLongFull <- dataLong(dataShort=dataSet, timeColumn=as.character(survModelFormula) [2], eventColumn=as.character(censModelFormula) [2])
#   }
#   
#   for(i in 1:length(trainIndices)) {
#     
#     # 0. Extract training, test data and responses
#     TrainSet <- dataSet [trainIndices [[i]], ]
#     if(length(trainIndices)!=1) {
#       TestSet <- dataSet [-trainIndices [[i]], ]
#     }
#     else {
#       TestSet <- TrainSet
#     }
#     
#     # 1. Convert training data to long format
#     if(!is.null(idColumn)) {
#       TrainLong <- dataLongTimeDep(dataSemiLong=TrainSet, timeColumn=as.character(survModelFormula) [2], eventColumn=as.character(censModelFormula) [2], idColumn=idColumn)
#     }
#     else {
#       TrainLong <- dataLong(dataShort=TrainSet, timeColumn=as.character(survModelFormula) [2], eventColumn=as.character(censModelFormula) [2])
#     }
#     
#     # 2. Convert response in training data to censoring variable
#     TrainLong <- dataCensoring(dataShort=TrainLong, respColumn="y", idColumn="obj")
#     
#     # 3. Convert test data to long format
#     if(!is.null(idColumn)) {
#       TestLong <- dataLongTimeDep(dataSemiLong=TestSet, timeColumn=as.character(survModelFormula) [2], eventColumn=as.character(censModelFormula) [2], idColumn=idColumn)
#     }
#     else {
#       TestLong <- dataLong (dataShort=TestSet, timeColumn=as.character(survModelFormula) [2], eventColumn=as.character(censModelFormula) [2])
#     }
#     
#     # 4. Fit censoring model on training data in long format
#     CensnewFormula <- update(censModelFormula, yCens ~ timeInt + .)
#     CensFit <- glm (formula=CensnewFormula, data=TrainLong, family=binomial(link=linkFunc), control=glm.control(maxit=2500))
#     
#     # 5. Estimate censoring survival curves inputs on test data
#     # Exclude cases with new factor levels in test data in long format
#     Check <- "error" %in% class(tryCatch(predict(CensFit, TestLong, type="response"), error= function (e) e))
#     if(Check) {
#       
#       # Which columns are factors in test data?
#       IndexFactor <- which(sapply(1:dim(TestLong)[2], function (x) is.factor(TestLong [, x]))==TRUE)
#       
#       # What are the levels of these factors?
#       TestLevelsFactor <- sapply(IndexFactor, function (x) levels(TestLong [, x]))
#       
#       # First column does not count (censoring process)
#       # What are the levels of the corresponding factors in the training data?
#       TrainLevelsFactor <- sapply(IndexFactor, function (x) levels(TrainLong [, x+1]))
#       
#       # Which levels of the test data exists in the training data factors?
#       InLevelsFactor <- lapply(1:length(TestLevelsFactor), function (x) which((TestLevelsFactor [[x]] %in% TrainLevelsFactor [[x]])==FALSE))
#       ExcludeRows <- lapply (1:length(IndexFactor), function (j) which(TestLong [, IndexFactor [j]] %in% TestLevelsFactor [[j]] [InLevelsFactor [[j]]]))
#       ExcludeRows <- do.call(c, ExcludeRows)
#       
#       # Convert Indices of left out test data to complete data set in long format
#       ExcludeRowsConv <- vector("integer", length(ExcludeRows))
#       for(j in 1:length(ExcludeRows)) {
#         I1 <- sapply(1:dim(TrainLongFull) [1], function(x) TrainLongFull [x, -1] == TestLong [ExcludeRows [j], -1])
#         ExcludeRowsConv [j] <- which(sapply(1:dim(I1) [2], function (x) all(I1 [, x]))==TRUE)
#       }
#       ExcludeRowsCensList [[i]] <- ExcludeRowsConv
#       
#       # Exclude rows of TestLong
#       TestLong <- TestLong [-ExcludeRows, ]
#     }
#     oneMinusLambdaList [[i]] <- 1 - predict(CensFit, TestLong, type="response")
#     
#     # 7. Estimate survival function inputs on test data
#     SurvnewFormula <- update(survModelFormula, y ~ timeInt + .)
#     SurvFit <- glm (formula=SurvnewFormula, data=TrainLong, family=binomial(link=linkFunc), control=glm.control(maxit=2500))
#     oneMinusLambdaSurvList [[i]] <- 1 - predict(SurvFit, TestLong, type="response")
#   }
#   
#   # 8. Estimate prediction error curve
#   # Estimate survival function of censoring process (complete data set)
#   oneMinuslambda <- do.call(c, oneMinusLambdaList)
#   # Merge excluded rows
#   ExcludeRowsCens <- do.call(c, ExcludeRowsCensList)
#   ExcludeRowsDataSet <- do.call(c, ExcludeRowsDataSetList)
#   if(!is.null(ExcludeRowsCens)) {
#     TrainLongFullExc <- TrainLongFull [-ExcludeRowsCens, ]
#   }
#   else {
#     TrainLongFullExc <- TrainLongFull
#   }
#   G <- aggregate(oneMinuslambda ~ obj, FUN = cumprod, data = TrainLongFullExc, simplify = FALSE)
#   if(!is.null(ExcludeRowsDataSet)) {
#     newEvent <- dataSet [-ExcludeRowsDataSet, as.character(censModelFormula) [2]]
#     newTime <- dataSet [-ExcludeRowsDataSet, as.character(survModelFormula) [2]]
#   }
#   else {
#     newEvent <- dataSet [, as.character(censModelFormula) [2]]
#     newTime <- dataSet [, as.character(survModelFormula) [2]]
#   }
#   n <- length(newEvent)
#   GT <- sapply(1:n, function(u){
#     if (newTime[u] > 1)
#       return(G[[2]] [u] [[1]] [newTime[u] - 1]) else
#         return(1) } )
#   
#   # Iterate over all time points
#   Lengthtimepoints <- length(timepoints)
#   weights <- vector("list", Lengthtimepoints)
#   results <- vector("numeric", Lengthtimepoints)
#   for(k in 1:Lengthtimepoints) {
#   
#   GTfixed <- sapply(1:n, function(u){
#   if (timepoints [k] > 1)
#     return(G[[2]] [u] [[1]] [timepoints [k] ]) else
#       return(1) } )
# 
#   # Estimate survival function
#   oneMinuslambdaSurv <- do.call(c, oneMinusLambdaSurvList)
#   S <- aggregate(oneMinuslambdaSurv ~ obj, FUN = cumprod, data = TrainLongFullExc, simplify = FALSE)
#   STfixed <- sapply(1:n, function(u){
#   if (timepoints [k] > 1)
#     return(S[[2]] [u] [[1]] [timepoints [k] ]) else
#       return(1) } )
#   
#   # Observed survival function must be evaluated on test data!
#   # Calculate observed survival function (complete data)
#   TrainLongFullExcSplit <- split(TrainLongFullExc$y, TrainLongFullExc$obj)
#   Sobs <- lapply(1:length(TrainLongFullExcSplit), ObsSurvFunc)
#   Sobs <- do.call(c, Sobs)
# 
#   # Exclude cases for which observed survival function is not available
#   IndexExclusion <- is.na(Sobs) | is.nan(Sobs) | is.infinite(Sobs)
#   if(!all(IndexExclusion==FALSE)) {
#     Sobs <- Sobs [!IndexExclusion]
#     STfixed <- STfixed [!IndexExclusion]
#     GTfixed <- GTfixed [!IndexExclusion]
#     GT <- GT [!IndexExclusion]
#     newEvent <- newEvent [!IndexExclusion]
#   }
#   
#   # Calculate prediction error curve given timepoints t
#   weights [[k]] <- WeightFunction ()
#   results [k] <- predErr ()
#   }
#   
#   # Combine outputs
#   RET <- list(Output=list(predErr = results, weights = weights), 
#               Input=list(timepoints=timepoints, dataSet=dataSet, trainIndices=trainIndices, survModelFormula=survModelFormula, censModelFormula=censModelFormula, linkFunc=linkFunc, idColumn=idColumn, Short=FALSE))
#   class(RET) <- "discSurvPredErrDisc"
#   return(RET)
# }

#' @rdname predErrCurve
#' @param x Object of class "discSurvPredErrDisc"("class discSurvPredErrDisc")
# #' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @keywords survival
#' @method print discSurvPredErrDisc
#' @export
print.discSurvPredErrDisc <- function (x, ...) {
  print(round(x$Output$predErr, 4))
}

#' @rdname predErrCurve
#' @param x Object of class "discSurvPredErrDisc"("class discSurvPredErrDisc")
#' @param \dots Additional arguments to S3 methods.
# #' @author Thomas Welchowski
#' @keywords survival
#' @method plot discSurvPredErrDisc
#' @export
plot.discSurvPredErrDisc <- function (x, ...) {
  
  plot(x = x$Input$timepoints, y = round(x$Output$predErr, 4), 
       xlab = "Time intervals", main = "Prediction error curve", ylab = "Value",
       type = "l", lty = 1, las = 1, ylim = c(0, 0.3), ...)
  abline(h = 0.25, lty = 2)

}

#####################################
# predErrCurve
# Short Version without CV but usable for arbitrary models



#' Prediction Error Curves
#' 
#' Estimates prediction error curves of arbitrary discrete survival prediction models. 
#' In prediction error curves the estimated and observed survival functions are
#' compared adjusted by weights at given timepoints.
#' 
#' The prediction error curves should be smaller than 0.25 for all time points,
#' because this is equivalent to a random assignment error.
#' 
#' @param timepoints Vector of the number of discrete time intervals ("integer vector").
#' @param estSurvList List of persons in the test data ("class list"). Each element contains a
#' estimated survival functions of all given time points ("numeric vector").
#' @param testTime Discrete survival times in the test data ("numeric vector").
#' @param testEvent Univariate event indicator in the test data ("binary vector").
#' @param trainTime Numeric vector of discrete survival times in the training
#' data ("numeric vector").
#' @param trainEvent Integer vector of univariate event indicator in the
#' training data("integer vector").
#' @return \itemize{ \item{List} List with objects: \itemize{ \item{Output} List
#' with two components \itemize{ \item{predErr} Numeric vector with estimated
#' prediction error values.  Names give the evaluation time point.
#' \item{weights} List of weights used in the estimation. Each list component
#' gives the weights of a person in the test data.  } \item{Input} A list of
#' given argument input values (saved for reference) } }
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link[mgcv]{gam}}
#' @references 
#' \insertRef{laanUniCensor}{discSurv} \cr\cr
#' \insertRef{gerdsConsisEst}{discSurv}
#' @keywords survival
#' @examples
#' 
#' # Example with cross validation and unemployment data 
#' library(Ecdat)
#' library(mgcv)
#' data(UnempDur)
#' summary(UnempDur$spell)
#' 
#' # Extract subset of data
#' set.seed(635)
#' IDsample <- sample(1:dim(UnempDur)[1], 100)
#' UnempDurSubset <- UnempDur [IDsample, ]
#' head(UnempDurSubset)
#' range(UnempDurSubset$spell)
#' 
#' # Generate training and test data
#' set.seed(7550)
#' TrainIndices <- sample (x = 1:dim(UnempDurSubset) [1], size = 75)
#' TrainUnempDur <- UnempDurSubset [TrainIndices, ]
#' TestUnempDur <- UnempDurSubset [-TrainIndices, ]
#' 
#' # Convert to long format
#' LongTrain <- dataLong(dataShort = TrainUnempDur, timeColumn = "spell", eventColumn = "censor1")
#' LongTest <- dataLong(dataShort = TestUnempDur, timeColumn = "spell", eventColumn = "censor1")
#' # Convert factor to numeric for smoothing
#' LongTrain$timeInt <- as.numeric(as.character(LongTrain$timeInt))
#' LongTest$timeInt <- as.numeric(as.character(LongTest$timeInt))
#' 
#' ######################################################################
#' # Estimate a generalized, additive model in discrete survival analysis
#' 
#' gamFit <- gam (formula = y ~ s(timeInt) + age + logwage, data = LongTrain, family = binomial())
#' 
#' # Estimate survival function of each person in the test data
#' oneMinusPredHaz <- 1 - predict(gamFit, newdata = LongTest, type = "response")
#' predSurv <- aggregate(oneMinusPredHaz ~ obj, data = LongTest, FUN = cumprod)
#' 
#' # Prediction error in first interval
#' tryPredErrDisc1 <- predErrCurve (timepoints = 1, 
#' estSurvList = predSurv [[2]], testTime = TestUnempDur$spell,
#' testEvent=TestUnempDur$censor1, trainTime = TrainUnempDur$spell,
#'  trainEvent=TrainUnempDur$censor1)
#' tryPredErrDisc1
#' 
#' # Prediction error of the 2. to 10. interval
#' tryPredErrDisc2 <- predErrCurve (timepoints = 2:10,
#' estSurvList = predSurv [[2]], testTime = TestUnempDur$spell,
#' testEvent = TestUnempDur$censor1, trainTime = TrainUnempDur$spell,
#' trainEvent = TrainUnempDur$censor1)
#' tryPredErrDisc2
#' plot(tryPredErrDisc2)
#' 
#' ########################################
#' # Fit a random discrete survival forest
#' 
#' library(ranger)
#' LongTrainRF <- LongTrain
#' LongTrainRF$y <- factor(LongTrainRF$y)
#' rfFit <- ranger(formula = y ~ timeInt + age + logwage, data = LongTrainRF,
#' probability = TRUE)
#' 
#' # Estimate survival function of each person in the test data
#' oneMinusPredHaz <- 1 - predict(rfFit, data = LongTest)$predictions[, 2]
#' predSurv <- aggregate(oneMinusPredHaz ~ obj, data = LongTest, FUN = cumprod)
#' 
#' # Prediction error in first interval
#' tryPredErrDisc1 <- predErrCurve (timepoints = 1, 
#' estSurvList = predSurv [[2]], testTime = TestUnempDur$spell,
#' testEvent = TestUnempDur$censor1, trainTime = TrainUnempDur$spell,
#'  trainEvent = TrainUnempDur$censor1)
#' tryPredErrDisc1
#' 
#' # Prediction error of the 2. to 10. interval
#' tryPredErrDisc2 <- predErrCurve (timepoints = 2:10,
#' estSurvList = predSurv [[2]], testTime = TestUnempDur$spell,
#' testEvent = TestUnempDur$censor1, trainTime = TrainUnempDur$spell,
#' trainEvent = TrainUnempDur$censor1)
#' tryPredErrDisc2
#' plot(tryPredErrDisc2)
#' 
#' @export predErrCurve
predErrCurve <- function (timepoints, estSurvList, testTime, 
                              testEvent, trainTime, trainEvent) {
  
  # Help functions
  WeightFunction <- function () {
    PartialSum1 <- testEventTemp * (1 - Sobs) / GT
    PartialSum2 <- Sobs / GTfixed
    return(PartialSum1 + PartialSum2)
  }
  predErr <- function () {
    sum(weights1[IncludeInd][finiteCheck] * 
          (estSurvInd[IncludeInd][finiteCheck] - 
             Sobs[IncludeInd][finiteCheck])^2) / sum(finiteCheck)
  }
  
  #####################
  # Execution
  
  # Expand training data in long format with censoring variable
  dataSetLong <- dataLong(dataShort = data.frame(
    trainTime = trainTime, trainEvent = trainEvent), 
    timeColumn = "trainTime", eventColumn = "trainEvent", timeAsFactor=TRUE)
  dataSetLongCens <- dataCensoring(dataShort = dataSetLong, timeColumn = "timeInt",
                                    shortFormat = FALSE)
  dataSetLongCens <- na.omit(dataSetLongCens)
  
  # Estimate glm with no covariates of censoring process
  glmCovariateFree <- glm(yCens ~ timeInt, data = dataSetLongCens, 
                          family = binomial(), control = glm.control(maxit = 2500))
  # Predict survival function of censoring process
  factorPrep <- factor(1:max(as.numeric(as.character(dataSetLongCens$timeInt))))
  GT_est <- cumprod(1 - predict(glmCovariateFree, 
                                newdata = data.frame(timeInt = factorPrep), 
                                type = "response"))
  
  # Loop over all time points
  predErrorValues <- vector("numeric", length(timepoints))
  StoreWeights <- vector("list", length(timepoints))
  for( k in 1:length(timepoints) ) {
    
    # Instable!
    # Restrict testTime and testEvent
    #   Check <- testTime < timepoints [k] & testEvent == 0
    #   testTimeTemp <- testTime [!Check]
    #   testEventTemp <- testEvent [!Check]
    #   if(length(testTimeTemp)==0) {
    #     return(0)
    #   }
    
    testTimeTemp <- testTime
    testEventTemp <- testEvent
    
    # Exclude observations, which were censored before t
    IncludeInd <- ifelse(testTimeTemp < timepoints [k] & 
                           testEventTemp == 0, FALSE, TRUE)
    
    # Estimate survival functions of censoring process
    GT <- c(1, GT_est)
    GT <- GT [testTimeTemp]
    GTfixed <- GT_est [timepoints [k] ]
    
    # Generate observed survival function
    Sobs <- ifelse(timepoints [k] < testTimeTemp, 1, 0)
    
    # Filter all predicted values: First order is timepoint and second layer is individuals
    estSurvInd <- sapply(1:length(estSurvList), 
                         function (x) estSurvList [[x]] [timepoints [k] ])
    # estSurvInd <- estSurvInd [!Check]
    
    # Estimate weights of each person in test data
    weights1 <- WeightFunction ()
    StoreWeights [[k]] <- weights1
    
    # Estimate prediction error
    finiteCheck <- is.finite(weights1[IncludeInd]) & 
      is.finite(estSurvInd[IncludeInd]) & 
      is.finite(Sobs[IncludeInd])
    predErrorValues [k] <- predErr ()
  }
  names(predErrorValues) <- paste("T=", timepoints, sep = "")
  # CheckRemove <- is.infinite(predErrorValues) | is.nan(predErrorValues) | is.na (predErrorValues)
  # CheckInclude <- is.finite(predErrorValues)
  # predErrorValues <- predErrorValues [CheckInclude]
  
  # Combine outputs
  RET <- list(Output=list(predErr = predErrorValues, weights = StoreWeights), 
              Input = list(timepoints = timepoints, estSurvList = estSurvList, 
                           testTime = testTime, testEvent = testEvent, 
                           trainTime = trainTime, trainEvent = trainEvent, Short = TRUE))
  class(RET) <- "discSurvPredErrDisc"
  return(RET)
}

#######################
# intPredErrCurve

# Description
# Estimates prediction error curves



#' Integrated Prediction Error Curves
#' 
#' Estimates the integrated prediction error curve based on a discrete survival model. 
#' This measure is time invariant and is based on the weighted quadratic differences of 
#' estimated and observed survival functions.
#' 
#' 
#' @param predErrObj Object generated by function \code{\link{predErrCurve}} ("class discSurvPredErrDisc").
#' @param tmax Gives the maximum time interval for which prediction errors are
#' calculated ("integer vector"). It must be smaller than the maximum observed time in the
#' training data of the object produced by function \code{\link{predErrCurve}}. 
#' Default value of NULL means, that all observed intervals are used.
#' @return Integrated prediction error per time interval ("numeric vector").
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link{predErrCurve}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv} \cr\cr
#' \insertRef{gneitingPropScore}{discSurv}
#' @keywords internal
#' @examples
#' 
#' #####################################################
#' # Example with cross validation and unemployment data 
#' 
#' library(Ecdat)
#' library(mgcv)
#' data(UnempDur)
#' summary(UnempDur$spell)
#' 
#' # Extract subset of data
#' set.seed(635)
#' IDsample <- sample(1:dim(UnempDur)[1], 100)
#' UnempDurSubset <- UnempDur [IDsample, ]
#' head(UnempDurSubset)
#' range(UnempDurSubset$spell)
#' 
#' # Generate training and test data
#' set.seed(7550)
#' TrainIndices <- sample (x = 1:dim(UnempDurSubset) [1], size = 75)
#' TrainUnempDur <- UnempDurSubset [TrainIndices, ]
#' TestUnempDur <- UnempDurSubset [-TrainIndices, ]
#' 
#' # Convert to long format
#' LongTrain <- dataLong(dataShort = TrainUnempDur, timeColumn = "spell", eventColumn = "censor1")
#' LongTest <- dataLong(dataShort = TestUnempDur, timeColumn = "spell", eventColumn = "censor1")
#' # Convert factor to numeric for smoothing
#' LongTrain$timeInt <- as.numeric(as.character(LongTrain$timeInt))
#' LongTest$timeInt <- as.numeric(as.character(LongTest$timeInt))
#' 
#' # Estimate, for example, a generalized, additive model in discrete survival analysis
#' gamFit <- gam (formula = y ~ s(timeInt) + age + logwage, data = LongTrain, family = binomial())
#' 
#' # Estimate survival function of each person in the test data
#' oneMinusPredHaz <- 1 - predict(gamFit, newdata = LongTest, type = "response")
#' predSurv <- aggregate(oneMinusPredHaz ~ obj, data = LongTest, FUN = cumprod)
#' 
#' # Prediction error in first and second interval
#' tryPredErrDisc1 <- predErrCurve (timepoints = 1, 
#' estSurvList = predSurv [[2]], newTime = TestUnempDur$spell,
#'  newEvent = TestUnempDur$censor1, trainTime = TrainUnempDur$spell,
#'  trainEvent = TrainUnempDur$censor1)
#' 
#' # Estimate integrated prediction error
#' tryintPredErrCurve <- intPredErrCurve (tryPredErrDisc1)
#' tryintPredErrCurve
#' 
#' # Example up to interval 3 (higher intervals are truncated)
#' tryintPredErrCurve2 <- intPredErrCurve (tryPredErrDisc1, tmax = 3)
#' tryintPredErrCurve2
#' 
#' ##########################
#' # Example with cancer data
#' 
#' library(survival)
#' head(cancer)
#' 
#' # Data preparation and convertion to 30 intervals
#' cancerPrep <- cancer
#' cancerPrep$status <- cancerPrep$status-1
#' intLim <- quantile(cancerPrep$time, prob = seq(0, 1, length.out = 30))
#' intLim [length(intLim)] <- intLim [length(intLim)] + 1
#' 
#' # Cut discrete time in smaller number of intervals
#' cancerPrep <- contToDisc(dataShort = cancerPrep, timeColumn = "time", intervalLimits = intLim)
#' 
#' # Generate training and test data
#' set.seed(753)
#' TrainIndices <- sample (x = 1:dim(cancerPrep) [1], size = dim(cancerPrep) [1] * 0.75)
#' TrainCancer <- cancerPrep [TrainIndices, ]
#' TestCancer <- cancerPrep [-TrainIndices, ]
#' TrainCancer$timeDisc <- as.numeric(as.character(TrainCancer$timeDisc))
#' TestCancer$timeDisc <- as.numeric(as.character(TestCancer$timeDisc))
#' 
#' # Convert to long format
#' LongTrain <- dataLong(dataShort = TrainCancer, timeColumn = "timeDisc", eventColumn = "status")
#' LongTest <- dataLong(dataShort = TestCancer, timeColumn = "timeDisc", eventColumn = "status")
#' # Convert factors
#' LongTrain$timeInt <- as.numeric(as.character(LongTrain$timeInt))
#' LongTest$timeInt <- as.numeric(as.character(LongTest$timeInt))
#' LongTrain$sex <- factor(LongTrain$sex)
#' LongTest$sex <- factor(LongTest$sex)
#' 
#' # Estimate, for example, a generalized, additive model in discrete survival analysis
#' gamFit <- gam (formula = y ~ s(timeInt) + s(age) + sex + ph.ecog, data = LongTrain, 
#' family = binomial())
#' summary(gamFit)
#' 
#' # Estimate survival function of each person in the test data
#' oneMinusPredHaz <- 1 - predict(gamFit, newdata = LongTest, type = "response")
#' predSurv <- aggregate(oneMinusPredHaz ~ obj, data = LongTest, FUN = cumprod)
#' 
#' # One survival curve is missing: Replace the missing values,
#' # with average value of other survival curves
#' predSurvTemp <- predSurv [[2]]
#' for(i in 1:length(predSurv [[2]])) {
#'   lenTemp <- length(predSurv [[2]] [[i]])
#'   if(lenTemp < 32) {
#'     predSurvTemp [[i]] <- c(predSurv [[2]] [[i]], rep(NA, 30 - lenTemp))
#'   }
#' }
#' # Calculate average survival curve
#' avgSurv <- rowMeans(do.call(cbind, predSurvTemp), na.rm=TRUE) [1:4]
#' # Insert into predicted survival curves
#' predSurvTemp <- predSurv [[2]]
#' for(i in 3:(length(predSurvTemp)+1)) {
#'   if(i == 3) {
#'     predSurvTemp [[i]] <- avgSurv
#'   }
#'   else {
#'     predSurvTemp [[i]] <- predSurv [[2]] [[i - 1]]
#'   }
#' }
#' # Check if the length of all survival curves is equal to the observed
#' # time intervals in test data
#' all(sapply(1:length(predSurvTemp), function (x) length(predSurvTemp [[x]]))==
#' as.numeric(as.character(TestCancer$timeDisc))) # TRUE
#' 
#' # Prediction error second interval
#' tryPredErrDisc1 <- predErrCurve (timepoints = 2, 
#' estSurvList = predSurvTemp, newTime = TestCancer$timeDisc,
#'  newEvent = TestCancer$status, trainTime = TrainCancer$timeDisc,
#'  trainEvent = TrainCancer$status)
#' 
#' # Calculate integrated prediction error
#' tryintPredErrCurve <- intPredErrCurve (tryPredErrDisc1)
#' tryintPredErrCurve
#' 
#' @noRd
intPredErrCurve <- function (predErrObj, tmax = NULL) {
  
  # Input check
  if(!(class(predErrObj) == "discSurvPredErrDisc")) {
    stop("Object *predErrObj* is not of class *discSurvPredErrDisc*! Please give an appropriate objecte type as input.")}
  
#  if(predErrObj$Input$Short==FALSE) {
#  
#   # Help function
#   predErrDiscTime <- function (t) {
#     predErrDisc (timepoints=t, 
#                  dataShort=predErrObj$Input$dataShort, trainIndices=predErrObj$Input$trainIndices, survModelFormula=predErrObj$Input$survModelFormula, censModelFormula=predErrObj$Input$censModelFormula, linkFunc=predErrObj$Input$linkFunc, idColumn=predErrObj$Input$idColumn)
#   }
# 
#   # Estimate marginal probabilities P(T=t)
#   MaxTime <- max(predErrObj$Input$dataShort [, as.character(predErrObj$Input$survModelFormula) [2] ])
#   DataSet <- predErrObj$Input$dataShort
#   TrainIndices <- predErrObj$Input$trainIndices
#   SurvModelFormula <- predErrObj$Input$survModelFormula
#   CensModelFormula <- predErrObj$Input$censModelFormula
#   LinkFunc <- predErrObj$Input$linkFunc
#   IdColumn <- predErrObj$Input$idColumn
# 
#   # Estimate survival function and marginal probabilities without covariates
#   if(!is.null(IdColumn)) {
#     TrainLongFull <- dataLongTimeDep(dataSemiLong=DataSet, timeColumn=as.character(SurvModelFormula) [2], eventColumn=as.character(CensModelFormula) [2], idColumn=IdColumn)
#   }
#   else {
#     TrainLongFull <- dataLong(dataShort=DataSet, timeColumn=as.character(SurvModelFormula) [2], eventColumn=as.character(CensModelFormula) [2])
#   }
#   MargFormula <- y ~ timeInt
#   MargFit <- glm (formula=MargFormula, data=TrainLongFull, family=binomial(link=LinkFunc), control=glm.control(maxit=2500))
#   PredMargData <- data.frame(timeInt=factor(min(TrainLongFull [, as.character(SurvModelFormula) [2] ]):max(TrainLongFull [, as.character(SurvModelFormula) [2] ])))
#   MargHaz <- as.numeric(predict(MargFit, PredMargData, type="response"))
#   # Marginal Probability P(T=t) = \lambda (T=t) \prod_{j=1}^{t-1} (1-\lambda (T=j))
#   MargProbs <- estMargProb(MargHaz)
# 
#   # Estimate prediction error curve over all time points
#   PredsErrorCurve <- predErrDiscTime (1:(MaxTime+1))$Output$predErr
#   IncludedIndices <- !(is.na(PredsErrorCurve) | is.nan (PredsErrorCurve) | is.infinite(PredsErrorCurve))
#   PredsErrorCurve <- PredsErrorCurve [IncludedIndices]
#   MargProbs <- as.numeric(MargProbs [IncludedIndices])
#   
#   # Output
#   Result <- sum(PredsErrorCurve * MargProbs)
#   return(c(IntPredErr=Result))
#   }
#  else {
    # Help function
    predErrDiscTime <- function (t) {
      predErrCurve (timepoints = t, estSurvList = EstSurvList, 
                        testTime = NewTime, testEvent = NewEvent, 
                        trainTime = TrainTime, trainEvent = TrainEvent)
    }
    
    # Estimate marginal probabilities P(T=t)
    MaxTime <- max(predErrObj$Input$trainTime)
    if(!is.null(tmax)) {
      if(tmax <= MaxTime) {
        MaxTime <- tmax
      }
      else {
        warning("Argument *tmax* is higher than the latest observed interval in training data. 
                Only prediction errors up to the latest observed interval time are given.")
      }
    }
    EstSurvList <- predErrObj$Input$estSurvList
    NewTime <- predErrObj$Input$newTime
    NewEvent <- predErrObj$Input$newEvent
    TrainTime <- predErrObj$Input$trainTime
    TrainEvent <- predErrObj$Input$trainEvent
    
    # Estimate survival function and marginal probabilities without covariates
    TrainLongFull <- dataLong(dataShort = data.frame(TrainTime = TrainTime, TrainEvent = TrainEvent), timeColumn = "TrainTime", 
                              eventColumn = "TrainEvent", timeAsFactor=TRUE)
    MargFormula <- y ~ timeInt
    MargFit <- glm (formula = MargFormula, data = TrainLongFull, family = binomial(), control = glm.control(maxit = 2500))
    PredMargData <- data.frame(timeInt = factor(1:MaxTime))
    MargHaz <- as.numeric(predict(MargFit, PredMargData, type = "response"))
    # Marginal Probability P(T=t) = \lambda (T=t) \prod_{j=1}^{t-1} (1-\lambda (T=j))
    MargProbs <- estMargProb(MargHaz)
    
    # Estimate prediction error curve over all time points
    PredsErrorCurve <- predErrDiscTime (1:MaxTime)$Output$predErr
    # IncludedIndices <- !(is.na(PredsErrorCurve) | is.nan (PredsErrorCurve) | is.infinite(PredsErrorCurve))
    IncludedIndices <- is.finite(PredsErrorCurve) & is.finite(MargProbs)
    PredsErrorCurve <- PredsErrorCurve [IncludedIndices]
    MargProbs <- as.numeric(MargProbs [IncludedIndices])
    
    # # Enlarge PredsErrorCurve with 0 to match length of MargProbs
    # if(length(PredsErrorCurve)!=length(MargProbs)) {
    #   PredsErrorCurve <- c(PredsErrorCurve, rep(0, length(MargProbs) - length(PredsErrorCurve)))
    # }
    
    # Output
    # Result <- sum(PredsErrorCurve [1:MaxTime] * MargProbs [1:MaxTime]) / sum(MargProbs [1:MaxTime])
    Result <- sum(PredsErrorCurve * MargProbs) / sum(MargProbs)
    return(c(IntPredErr = Result))
#  }
}

######################
# martingaleResid

# Description
# Calculates the martingale residuals


#' @title Martingale Residuals
#' 
#' @description Estimates the martingale residuals of discrete survival model.
#' 
#' @param hazards Predicted hazards from a discrete survival model ("numeric vector").
#' @param dataSetLong Data set in long format ("class data.frame").
#' @return Martingale residuals for each observation in long format ("numeric vector").
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link[stats]{glm}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv} \cr\cr
#' \insertRef{therneauMart}{discSurv}
#' @keywords internal
#' @examples
#' 
#' # Example with cross validation and unemployment data 
#' library(Ecdat)
#' data(UnempDur)
#' summary(UnempDur$spell)
#' 
#' # Extract subset of data
#' set.seed(635)
#' IDsample <- sample(1:dim(UnempDur)[1], 100)
#' UnempDurSubset <- UnempDur [IDsample, ]
#' 
#' # Conversion to long format
#' UnempDurSubsetLong <- dataLong(dataShort = UnempDurSubset,
#' timeColumn = "spell", eventColumn = "censor1")
#' 
#' # Estimate discrete survival continuation ratio model
#' contModel <- glm(y ~ timeInt + age + logwage, data = UnempDurSubsetLong,
#' family = binomial(link = "logit"))
#' 
#' # Fit hazards to the data set in long format
#' hazPreds <- predict(contModel, type = "response")
#' 
#' # Calculate martingale residuals for the unemployment data subset
#' MartResid <- martingaleResid (hazards = hazPreds, dataSetLong = UnempDurSubsetLong)
#' MartResid
#' sum(MartResid)
#' 
#' # Plot martingale residuals vs each covariate in the event interval
#' # Dotted line represents the loess estimate
#' plot(MartResid, covariates = c("age", "logwage"), dataSetLong = UnempDurSubsetLong)
#' 
#' @export martingaleResid
martingaleResid <- function(hazards, dataSetLong) {
  
  # Input checks
  if(!is.data.frame(dataSetLong)) {stop("Argument *dataSetLong* is not in the correct format! Please specify as data.frame object.")}

  # Calculate residuals
  splithazPreds <- split(hazards, dataSetLong$obj)
  splitY <- split(dataSetLong$y, dataSetLong$obj)
  martResid <- sapply(1:length(splitY), function (x) sum(splitY [[x]] - splithazPreds [[x]]))
  class(martResid) <- "discSurvMartingaleResid"
  return(martResid)
}

#' @rdname martingaleResid
#' @param x Object of class "discSurvMartingaleResid"("class discSurvMartingaleResid")
#' @param covariates Names of covariates to plot ("character vector").
#' @param dataSetLong Data in long format ("class data.frame").
#' @param \dots Additional arguments to the plot function
# #' @author Thomas Welchowski
#' @details Gives a different plot of each marginal covariate against the martingale
#' residuals. Additionally a nonparametric \code{\link{loess}} estimation is
#' done.
#' @keywords internal
#' @method plot discSurvMartingaleResid
#' @export
plot.discSurvMartingaleResid <- function (x, covariates, dataSetLong, ...) {

  splitDataSetLong <- split(dataSetLong, dataSetLong$obj)
  tailSplitDataSetLong <- lapply(splitDataSetLong, function (j) tail(j, 1))
  tailSplitDataSetLong <- do.call(rbind, tailSplitDataSetLong)
  
  # Covariates of interest
  lengthCovariates <- length(covariates)

  # Create plots of covariates
  for( i in 1:lengthCovariates ) {
    tempData <- data.frame(x = tailSplitDataSetLong [, covariates [i] ], y = c(x))
    tempData <- tempData [order(tempData$x), ]
    plot(x = tempData$x, y = tempData$y, las = 1, xlab = covariates [i], 
         ylab = "Martingale Residuals", ...)
    if(is.numeric(tempData$x)) {
      loessPred <- predict(loess(formula = y ~ x, data = tempData))
      lines(x = tempData$x, y = loessPred)
    }
    abline(h = 0, lty = 2)
  }
}

########################
# devResid

# Description
# Computes the root of the squared deviance residual



#' @title Deviance Residuals
#' 
#' @description Computes the root of the deviance residuals for evaluation of performance in
#' discrete survival analysis. A generalized, linear model is used for
#' prediction.
#' 
#' 
#' @param dataSet Original data in short format ("class data.frame").
#' @param survModelFormula Gives the specified relationship of discrete
#' response and covariates ("class formula"). The formula is designed, that the intercepts for
#' the time dependent base line hazards are always included. Therefore only
#' covariates should be given in this formula. This argument is required to be
#' of class "formula".
#' @param eventColumn Gives the column name of the event indicator (1=observed,
#' 0=censored) ("character vector").
#' @param linkFunc Specifies the desired link function in use of generalized,
#' linear models ("character vector").
#' @param idColumn Gives the column name of the identification number of each
#' person ("character vector"). Default value of NULL means, that
#' each row equals one person (no repeated measurements).
#' @return \itemize{ \item{Output} List with objects: \itemize{ \item{DevResid}
#' Square root of deviance residuals as numeric vector.  \item{GlmFit} Fit
#' object of class (generalized, linear model used in the calculations) }
#' \item{Input} A list of given argument input values (saved for reference) }
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link{adjDevResidGlm}}, \code{\link{glm}}, 
#' \code{\link{predErrCurve}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv} \cr\cr
#' \insertRef{tutzRegCat}{discSurv}
#' @keywords internal
#' @examples
#' 
#' library(Ecdat)
#' # Example with cross validation and unemployment data 
#' data(UnempDur)
#' summary(UnempDur$spell)
#' 
#' # Extract subset of data
#' set.seed(635)
#' IDsample <- sample(1:dim(UnempDur)[1], 100)
#' UnempDurSubset <- UnempDur [IDsample, ]
#' 
#' # Calculate deviance residuals for the unemployment data subset
#' devianceResiduals <- devResid (dataSet = UnempDurSubset, survModelFormula = spell ~ age + logwage, 
#' eventColumn = "censor1", linkFunc = "logit", idColumn = NULL)
#' devianceResiduals
#' 
#' @noRd
devResid <- function (dataSet, survModelFormula, eventColumn, linkFunc = "logit", idColumn = NULL) {
  
  # Input checks
  if(!is.data.frame(dataSet)) {stop("Argument *dataSet* is not in the correct format! Please specify as data.frame object.")}
  if(!("formula" %in% class(survModelFormula))) {stop("*survModelFormula* is not of class formula! Please specify a valid formula, e. g. y ~ x + z.")}
  if(!any(names(dataSet) == eventColumn)) {stop("Argument *eventColumn* is not available in *dataSet*! Please specify the correct column name of the event indicator.")}
  if(!(any(names(dataSet) == idColumn) | is.null(idColumn))) {stop("Argument *idColumn* is not available in *dataSet*! Please specify the correct column name of the identification numbers of persons.")}
  
  # Help function
  SquDevResid <- function (x) {-2 * sum(splitY [[x]] * log(splitHazards [[x]]) + (1 - splitY [[x]]) * log(1 - splitHazards [[x]] ))}
  
  # Convert to long format
  if(!is.null(idColumn)) {
    dataSetLong <- dataLongTimeDep(dataSemiLong = dataSet, timeColumn = as.character(survModelFormula) [2], eventColumn = eventColumn, idColumn = idColumn)
  }
  else {
    dataSetLong <- dataLong(dataShort = dataSet, timeColumn = as.character(survModelFormula) [2], eventColumn = eventColumn)
  }
  
  # Fit generalized, linear model
  NewFormula <- update(survModelFormula, y ~ timeInt + .)
  glmFit <- glm(formula = NewFormula, data = dataSetLong, family = binomial(link = linkFunc), control = glm.control(maxit = 2500))
  hazards <- predict(glmFit, type = "response")
  
  # Calculate residuals
  splitHazards <- split(hazards, dataSetLong$obj)
  splitY <- split(dataSetLong$y, dataSetLong$obj)
  Residuals <- sapply(1:length(splitY), SquDevResid)
  Output <- list(Output = list(DevResid = sqrt(Residuals), GlmFit = glmFit), 
                 Input = list(dataSet = dataSet, survModelFormula = survModelFormula, eventColumn = eventColumn, linkFunc = linkFunc, idColumn = idColumn))
  class(Output) <- "discSurvDevResid"
  return(Output)
}

#' @rdname devResid
#' @param x Object of class "discSurvDevResid"("class discSurvDevResis")
#' @param \dots Additional arguments to S3 methods.
# #' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @keywords internal
#' @method print discSurvDevResid
#' @noRd
print.discSurvDevResid <- function (x, ...) {
  print(round(x$Output$DevResid, 4))
}

########################
# adjDevResidGlm

# Description
# Calculates the adjusted deviance residuals. Should be normal distributed, in the case of a well fitting model.

#' @name adjDevResidGlm
#' @title Adjusted Deviance Residuals based on generalized linear models
#' 
#' @description Calculates the adjusted deviance residuals. This function only supports
#' generalized, linear models. The adjusted deviance residuals should be
#' approximately normal distributed, in the case of a well fitting model.
#' 
#' @param dataShort Original data in short format ("class data.frame").
#' @param survModelFormula Gives the specified relationship of discrete
#' response and covariates("class formula"). The formula is designed, that the intercepts for
#' the time dependent base line hazards are always included. Therefore only
#' covariates should be given in this formula. This argument is required to be
#' of class "formula".
#' @param eventColumn Gives the column name of the event indicator (1=observed,
#' 0=censored) ("character vector").
#' @param linkFunc Specifies the desired link function in use of generalized,
#' linear models("character vector").
#' @param idColumn Gives the column name of the identification number of each
#' person ("character vector"). Default value of NULL means, that
#' each row equals one person (no repeated measurements).
#' @param \dots Additional arguments to S3 methods.
#' @return \itemize{ \item{Output} List with objects: \itemize{
#' \item{AdjDevResid} Adjusted deviance residuals as numeric vector.
#' \item{GlmFit} Fit object of class (generalized, linear model used in the
#' calculations) } \item{Input} A list of given argument input values (saved for
#' reference) }
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link{devResid}}, \code{\link{glm}}, \code{\link{predErrCurve}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv} \cr\cr
#' \insertRef{tutzRegCat}{discSurv}
#' @keywords internal
#' @examples
#' 
#' # Example with cross validation and unemployment data 
#' library(Ecdat)
#' data(UnempDur)
#' summary(UnempDur$spell)
#' 
#' # Extract subset of data
#' set.seed(635)
#' IDsample <- sample(1:dim(UnempDur)[1], 100)
#' UnempDurSubset <- UnempDur [IDsample, ]
#' 
#' # Calculate adjusted deviance residuals for the unemployment data subset
#' adjDevianceResiduals <- adjDevResidGlm(dataShort = UnempDurSubset, 
#' survModelFormula = spell ~ age + logwage, 
#' eventColumn = "censor1", linkFunc = "logit", idColumn = NULL)
#' adjDevianceResiduals
#' 
#' # Exclude outliers
#' adjDevResidWoOut <- adjDevianceResiduals$Output$AdjDevResid [
#' adjDevianceResiduals$Output$AdjDevResid < 
#' quantile(adjDevianceResiduals$Output$AdjDevResid, prob=0.9) & 
#' adjDevianceResiduals$Output$AdjDevResid > 
#' quantile(adjDevianceResiduals$Output$AdjDevResid, prob=0.1)]
#' 
#' # Compare nonparametric density estimate of adjusted deviance residuals 
#' # with adapted normal distribution
#' plot(density(adjDevResidWoOut), xlim = c(-10,10), 
#' main = "Density comparison: Normal vs nonparametric estimate", lwd = 2, las = 1, col = "red")
#' dnorm1 <- function (x) {dnorm (x, mean = mean(adjDevResidWoOut), sd = sd(adjDevResidWoOut))}
#' curve (dnorm1, n = 500, add = TRUE, col = "black", lwd = 2)
#' legend("topright", legend = c("Normal", "Nonpar"), lty = 1, lwd = 2, col = c("black", "red"))
#' 
#' @noRd
adjDevResidGlm <- function (dataShort, survModelFormula, eventColumn, linkFunc = "logit", idColumn = NULL) {
  # Input checks
  if(!is.data.frame(dataShort)) {stop("Argument *dataShort* is not in the correct format! Please specify as data.frame object.")}
  if(!("formula" %in% class(survModelFormula))) {stop("*survModelFormula* is not of class formula! Please specify a valid formula, e. g. y ~ x + z.")}
  if(!any(names(dataShort) == eventColumn)) {stop("Argument *eventColumn* is not available in *dataShort*! Please specify the correct column name of the event indicator.")}
  if(!(any(names(dataShort) == idColumn) | is.null(idColumn))) {stop("Argument *idColumn* is not available in *dataShort*! Please specify the correct column name of the identification numbers of persons.")}
  
  # Help function
  AdjDevResid <- function (x) { 
    LogTerm1 <- ifelse(splitY [[x]] == 1, -log(splitHazards [[x]]), 0)
    LogTerm2 <- ifelse(splitY [[x]] == 0, -log (1 - splitHazards [[x]]), 0)
    FirstPartialSum <- sum(sign(splitY [[x]] - splitHazards [[x]]) * (sqrt(splitY [[x]] * LogTerm1 + (1 - splitY [[x]]) * LogTerm2)))
    SecondPartialSum <- sum( (1 - 2 * splitHazards [[x]]) / sqrt (splitHazards [[x]] * (1 - splitHazards [[x]]) * 36) )
    return(FirstPartialSum + SecondPartialSum)
  }
  
  # Convert to long format
  if(!is.null(idColumn)) {
    dataShortLong <- dataLongTimeDep(dataSemiLong = dataShort, timeColumn = as.character(survModelFormula) [2], eventColumn = eventColumn, idColumn = idColumn)
  }
  else {
    dataShortLong <- dataLong(dataShort = dataShort, timeColumn = as.character(survModelFormula) [2], eventColumn = eventColumn)
  }
  
  # Fit generalized, linear model
  NewFormula <- update(survModelFormula, y ~ timeInt + .)
  glmFit <- glm(formula = NewFormula, data = dataShortLong, family = binomial(link = linkFunc))
  hazards <- predict(glmFit, type = "response")
  
  # Calculate residuals
  splitHazards <- split(hazards, dataShortLong$obj)
  splitY <- split(dataShortLong$y, dataShortLong$obj)
  Residuals <- sapply(1:length(splitY), AdjDevResid)
  Output <- list(Output = list(AdjDevResid = Residuals, GlmFit = glmFit), 
                 Input = list(dataShort = dataShort, survModelFormula = survModelFormula, eventColumn = eventColumn, linkFunc = linkFunc, idColumn = idColumn))
  class(Output) <- "discSurvAdjDevResid"
  return(Output)
}

#' @rdname adjDevResidGlm
#' @param x Object of class "discSurvAdjDevResid"("class discSurvAdjDevResid")
# #' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @keywords internal
#' @method print discSurvAdjDevResid
#' @noRd
print.discSurvAdjDevResid <- function (x, ...) {
  print(round(x$Output$AdjDevResid, 4))
}

#' @rdname adjDevResidGlm
#' @param x Object of class "discSurvAdjDevResid"("class discSurvAdjDevResid")
# #' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @details Is called implicitly by using standard plot function on an object of class
#' "discSurvAdjDevResid". It plots a qqplot against the normal distribution. If
#' the model fits the data well, it should be approximately normal distributed.
#' @keywords internal
#' @method plot discSurvAdjDevResid
#' @noRd
plot.discSurvAdjDevResid <- function (x, ...) {
  qqnorm (y = x$Output$AdjDevResid, las = 1, ...)
  qqline(y = x$Output$AdjDevResid, ...)
}

########################
# adjDevResid

# Description
# Calculates the adjusted deviance residuals for arbitrary hazard prediction models. Should be normal distributed, in the case of a well fitting model.

#' Adjusted Deviance Residuals in short format
#' 
#' Calculates the adjusted deviance residuals for arbitrary prediction models.
#' The adjusted deviance residuals should be approximately normal distributed,
#' in the case of a well fitting model.
#' 
#' 
#' @param dataLong Data set in long format ("class data.frame").
#' @param hazards Estimated discrete hazards of the data in long format("numeric vector"). Hazard
#' rates are probabilities and therefore restricted to the interval [0, 1].
#' @return \itemize{ \item{Output} List with objects: \itemize{
#' \item{AdjDevResid} Adjusted deviance residuals as numeric vector }
#' \item{Input} A list of given argument input values (saved for reference) }
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link{intPredErr}}, \code{\link{predErrCurve}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv} \cr\cr
#' \insertRef{tutzRegCat}{discSurv}
#' @keywords survival
#' @examples
#' 
#' library(survival)
#' 
#' # Transform data to long format
#' heart[, "stop"] <- ceiling(heart[, "stop"])
#' set.seed(0)
#' Indizes <- sample(unique(heart$id), 25)
#' randSample <- heart[unlist(sapply(1:length(Indizes), 
#' function(x) which(heart$id == Indizes[x]))),]
#' heartLong <- dataLongTimeDep(dataSemiLong = randSample, 
#' timeColumn = "stop", eventColumn = "event", idColumn = "id", timeAsFactor = FALSE)
#' 
#' # Fit a generalized, additive model and predict discrete hazards on data in long format
#' library(mgcv)
#' gamFit <- gam(y ~ timeInt + surgery + transplant + s(age), data = heartLong, family = "binomial")
#' hazPreds <- predict(gamFit, type = "response")
#' 
#' # Calculate adjusted deviance residuals
#' devResiduals <- adjDevResid(dataLong = heartLong, hazards = hazPreds)$Output$AdjDevResid
#' devResiduals
#' 
#' @export adjDevResid
adjDevResid <- function(dataLong, hazards) {
  # Input checks
  if(!is.data.frame(dataLong)) {stop("Argument *dataLong* is not in the correct format! Please specify as data.frame object.")}
  if(!all(hazards >= 0 & hazards <= 1)) {stop("Argument *hazards* must contain probabilities in the closed interval from zero to one. Please verify that *hazards* are estimated discrete hazards")}
  if(!(dim(dataLong)[1] == length(hazards))) {stop("The length of argument *hazards* must match the number of observations")}
  
  # Help function
  AdjDevResid <- function (x) { 
    LogTerm1 <- ifelse(splitY [[x]] == 1, -log(splitHazards [[x]]), 0)
    LogTerm2 <- ifelse(splitY [[x]] == 0, -log (1 - splitHazards [[x]]), 0)
    FirstPartialSum <- sum(sign(splitY [[x]] - splitHazards [[x]]) * (sqrt(splitY [[x]] * LogTerm1 + (1 - splitY [[x]]) * LogTerm2)))
    SecondPartialSum <- sum( (1 - 2*splitHazards [[x]]) / sqrt (splitHazards [[x]] * (1 - splitHazards [[x]]) * 36) )
    return(FirstPartialSum + SecondPartialSum)
  }

  # Calculate residuals
  splitHazards <- split(hazards, dataLong$obj)
  splitY <- split(dataLong$y, dataLong$obj)
  Residuals <- sapply(1:length(splitY), AdjDevResid)
  Output <- list(Output = list(AdjDevResid = Residuals), 
                 Input = list(dataLong = dataLong, hazards = hazards))
  class(Output) <- "discSurvAdjDevResid"
  return(Output)
}

########################
# devResid

# Description
# Computes the root of the squared deviance residual



#' Deviance Residuals
#' 
#' Computes the root of the deviance residuals for evaluation of performance in
#' discrete survival analysis.
#' 
#' 
#' @param dataLong Original data in long format ("class data.frame").
#' The correct format can be specified with data preparation, see e. g.
#' \code{\link{dataLong}}.
#' @param hazards Estimated discrete hazards of the data in long format("numeric vector"). Discrete
#' discrete hazards are probabilities and therefore restricted to the interval [0,
#' 1].
#' @return \itemize{ \item{Output} List with objects: \itemize{ \item{DevResid}
#' Square root of deviance residuals as numeric vector.  } \item{Input} A list
#' of given argument input values (saved for reference) }
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link{adjDevResid}}, \code{\link{predErrCurve}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv} \cr\cr
#' \insertRef{tutzRegCat}{discSurv}
#' @keywords survival
#' @examples
#' 
#' library(survival)
#' 
#' # Transform data to long format
#' heart[, "stop"] <- ceiling(heart[, "stop"])
#' set.seed(0)
#' Indizes <- sample(unique(heart$id), 25)
#' randSample <- heart[unlist(sapply(1:length(Indizes), 
#' function(x) which(heart$id == Indizes[x]))),]
#' heartLong <- dataLongTimeDep(dataSemiLong = randSample, 
#' timeColumn = "stop", eventColumn = "event", idColumn = "id", timeAsFactor = FALSE)
#' 
#' # Fit a generalized, additive model and predict discrete hazards on data in long format
#' library(mgcv)
#' gamFit <- gam(y ~ timeInt + surgery + transplant + s(age), data = heartLong, family = "binomial")
#' hazPreds <- predict(gamFit, type = "response")
#' 
#' # Calculate the deviance residuals
#' devResiduals <- devResid (dataLong = heartLong, hazards = hazPreds)$Output$DevResid
#' 
#' # Compare with estimated normal distribution
#' plot(density(devResiduals), 
#' main = "Empirical density vs estimated normal distribution", 
#' las = 1, ylim = c(0, 0.5))
#' tempFunc <- function (x) dnorm(x, mean = mean(devResiduals), sd = sd(devResiduals))
#' curve(tempFunc, xlim = c(-10, 10), add = TRUE, col = "red")
#' # The empirical density seems like a mixture distribution, 
#' # but is not too far off in with values greater than 3 and less than 1
#' 
#' @export devResid
devResid <- function(dataLong, hazards) {
  
  # Input checks
  if(!is.data.frame(dataLong)) {stop("Argument *dataLong* is not in the correct format! Please specify as data.frame object.")}
  if(!all(hazards >= 0 & hazards<=1)) {stop("Argument *hazards* must contain probabilities in the closed interval from zero to one. Please verify that *hazards* are estimated discrete hazards")}
  if(!(dim(dataLong)[1] == length(hazards))) {stop("The length of argument *hazards* must match the number of observations")}
  
  # Help function
  SquDevResid <- function (x) {-2 * sum(splitY [[x]] * log(splitHazards [[x]]) + (1 - splitY [[x]]) * log(1 - splitHazards [[x]] ))}

  # Calculate residuals
  splitHazards <- split(hazards, dataLong$obj)
  splitY <- split(dataLong$y, dataLong$obj)
  Residuals <- sapply(1:length(splitY), SquDevResid)
  Output <- list(Output = list(DevResid = sqrt(Residuals)), 
                 Input = list(dataLong = dataLong, hazards = hazards))
  class(Output) <- "discSurvDevResid"
  return(Output)
}

##################
# C-Index

#' Concordance index
#' 
#' Calculates the concordance index for discrete survival models, which does not depend on time. 
#' This is the probability that, for a pair of randomly chosen comparable samples, the sample with the higher risk prediction will
#' experience an event before the other sample or belongs to a higher binary class.
#' 
#' 
#' @param marker Gives the predicted values of the linear predictor of a
#' regression model ("numeric vector"). May also be on the response scale.
#' @param testTime New time intervals in the test data ("integer vector").
#' @param testEvent Event indicators in the test data ("binary vector").
#' @param trainTime Time intervals in the training data ("integer vector").
#' @param trainEvent Event indicators in the training data ("binary vector").
#' @return Value of discrete concordance index between zero and one ("numeric
#' vector").
#' @note It is assumed that all time points up to the last observed interval
#' [a_{q-1}, a_q) are available.
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @references 
#' \insertRef{schmidDiscMeasure}{discSurv} \cr\cr
#' \insertRef{unoEvalPred}{discSurv} \cr\cr
#' \insertRef{heagertySurvROC}{discSurv}
#' @keywords survival
#' @examples
#' 
#' ##################################################
#' # Example with unemployment data and prior fitting
#' 
#' library(Ecdat)
#' library(caret)
#' library(mgcv)
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
#' UnempDurSubsetTrainLong <- dataLong(dataShort = UnempDurSubsetTrain, 
#' timeColumn = "spell", eventColumn = "censor1")
#' 
#' # Estimate gam with smooth baseline
#' gamFit <- gam(formula = y ~ s(I(as.numeric(as.character(timeInt)))) + 
#' s(age) + s(logwage), data = UnempDurSubsetTrainLong, family = binomial())
#' gamFitPreds <- predict(gamFit, newdata = cbind(UnempDurSubsetTest, 
#' timeInt = UnempDurSubsetTest$spell))
#' 
#' # Evaluate C-Index based on short data format
#' cIndex(marker = gamFitPreds, 
#' testTime = UnempDurSubsetTest$spell, 
#' testEvent = UnempDurSubsetTest$censor1, 
#' trainTime = UnempDurSubsetTrain$spell, 
#' trainEvent = UnempDurSubsetTrain$censor1)
#' 
#' #####################################
#' # Example National Wilm's Tumor Study
#' 
#' library(survival)
#' head(nwtco)
#' summary(nwtco$rel)
#' 
#' # Select subset
#' set.seed(-375)
#' Indices <- sample(1:dim(nwtco)[1], 500)
#' nwtcoSub <- nwtco [Indices, ]
#' 
#' # Convert time range to 30 intervals
#' intLim <- quantile(nwtcoSub$edrel, prob = seq(0, 1, length.out = 30))
#' intLim [length(intLim)] <- intLim [length(intLim)] + 1
#' nwtcoSubTemp <- contToDisc(dataShort = nwtcoSub, timeColumn = "edrel", intervalLimits = intLim)
#' nwtcoSubTemp$instit <- factor(nwtcoSubTemp$instit)
#' nwtcoSubTemp$histol <- factor(nwtcoSubTemp$histol)
#' nwtcoSubTemp$stage <- factor(nwtcoSubTemp$stage)
#' 
#' # Split in training and test sample
#' set.seed(-570)
#' TrainingSample <- sample(1:dim(nwtcoSubTemp)[1], round(dim(nwtcoSubTemp)[1]*0.75))
#' nwtcoSubTempTrain <- nwtcoSubTemp [TrainingSample, ]
#' nwtcoSubTempTest <- nwtcoSubTemp [-TrainingSample, ]
#' 
#' # Convert to long format
#' nwtcoSubTempTrainLong <- dataLong(dataShort = nwtcoSubTempTrain, 
#' timeColumn = "timeDisc", eventColumn = "rel", timeAsFactor=TRUE)
#' 
#' # Estimate glm
#' inputFormula <- y ~ timeInt + histol + instit + stage
#' glmFit <- glm(formula = inputFormula, data = nwtcoSubTempTrainLong, family = binomial())
#' linPreds <- predict(glmFit, newdata = cbind(nwtcoSubTempTest, 
#' timeInt = factor(nwtcoSubTempTest$timeDisc, levels=levels(nwtcoSubTempTrainLong$timeInt))))
#' 
#' # Evaluate C-Index based on short data format
#' cIndex(marker = linPreds, 
#' testTime = as.numeric(as.character(nwtcoSubTempTest$timeDisc)), 
#' testEvent = nwtcoSubTempTest$rel, 
#' trainTime = as.numeric(as.character(nwtcoSubTempTrain$timeDisc)), 
#' trainEvent = nwtcoSubTempTrain$rel) 
#' 
#' 
#' @export cIndex
cIndex <- function(marker, testTime, testEvent, trainTime, trainEvent){
  tprFit <- tprUnoShort(timepoint = 1, marker, testTime,
                         testEvent, trainTime, trainEvent)
  fprFit <- fprUnoShort(timepoint = 1, marker, testTime)
  aucFit <- aucUno(tprFit, fprFit)
  CFit <- unname(concorIndex(aucFit)$Output)
  return(CFit)
}

###########################################
# Integrated prediction error curve

# Integrated prediction error curves 


#' Integrated prediction error
#' 
#' Computes the integrated prediction error curve for discrete survival models.
#' 
#' 
#' @param hazards Predicted discrete hazards in the test data ("numeric vector").
#' @param testTime Discrete time intervals in short format of the test set
#' ("integer vector").
#' @param testEvent Events in short format in the test set ("binary vector").
#' @param trainTime Discrete time intervals in short format of the
#' training data set ("integer vector").
#' @param trainEvent Events in short format in the training set ("binary
#' vector").
#' @param testDataLong Test data in long format("class data.frame"). The discrete survival function is
#' calculated based on the predicted hazards. It is assumed that the data was
#' preprocessed with a function with prefix "dataLong", see e. g.
#' \code{\link{dataLong}}, \code{\link{dataLongTimeDep}}.
#' @param tmax Gives the maximum time interval for which prediction errors are
#' calculated ("integer vector"). It must be smaller than the maximum observed time in the
#' training data of the object produced by function. The default setting NULL means, that all observed intervals are used.
#' @return Integrated prediction error ("numeric vector").
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link{predErrCurve}}, \code{\link[stats]{aggregate}}
#' @references
#' \insertRef{tutzModelDisc}{discSurv} \cr\cr
#' \insertRef{gneitingPropScore}{discSurv}
#' @keywords survival
#' @examples
#' 
#' ##########################
#' # Example with cancer data
#' 
#' library(survival)
#' head(cancer)
#' 
#' # Data preparation and convertion to 30 intervals
#' cancerPrep <- cancer
#' cancerPrep$status <- cancerPrep$status-1
#' intLim <- quantile(cancerPrep$time, prob = seq(0, 1, length.out = 30))
#' intLim [length(intLim)] <- intLim [length(intLim)] + 1
#' 
#' # Cut discrete time in smaller number of intervals
#' cancerPrep <- contToDisc(dataShort = cancerPrep, timeColumn = "time", intervalLimits = intLim)
#' 
#' # Generate training and test data
#' set.seed(753)
#' TrainIndices <- sample (x = 1:dim(cancerPrep) [1], size = dim(cancerPrep) [1] * 0.75)
#' TrainCancer <- cancerPrep [TrainIndices, ]
#' TestCancer <- cancerPrep [-TrainIndices, ]
#' TrainCancer$timeDisc <- as.numeric(as.character(TrainCancer$timeDisc))
#' TestCancer$timeDisc <- as.numeric(as.character(TestCancer$timeDisc))
#' 
#' # Convert to long format
#' LongTrain <- dataLong(dataShort = TrainCancer, timeColumn = "timeDisc", eventColumn = "status",
#' timeAsFactor=FALSE)
#' LongTest <- dataLong(dataShort = TestCancer, timeColumn = "timeDisc", eventColumn = "status",
#' timeAsFactor=FALSE)
#' # Convert factors
#' LongTrain$timeInt <- as.numeric(as.character(LongTrain$timeInt))
#' LongTest$timeInt <- as.numeric(as.character(LongTest$timeInt))
#' LongTrain$sex <- factor(LongTrain$sex)
#' LongTest$sex <- factor(LongTest$sex)
#' 
#' # Estimate, for example, a generalized, additive model in discrete survival analysis
#' library(mgcv)
#' gamFit <- gam (formula = y ~ s(timeInt) + s(age) + sex + ph.ecog, data = LongTrain, 
#' family = binomial())
#' summary(gamFit)
#' 
#' # 1. Specification of predicted discrete hazards
#' # Estimate survival function of each person in the test data
#' testPredHaz <- predict(gamFit, newdata = LongTest, type = "response")
#' 
#' # 2. Calculate integrated prediction error
#' intPredErr(hazards = testPredHaz, 
#' testTime = TestCancer$timeDisc, testEvent = TestCancer$status, 
#' trainTime = TrainCancer$timeDisc, trainEvent = TrainCancer$status, 
#' testDataLong = LongTest)
#' 
#' @export intPredErr
intPredErr <- function(hazards, testTime, 
                           testEvent, trainTime, trainEvent, 
                           testDataLong, tmax = NULL){
  
    # Convert to 1-preds
    oneMinusPredHaz <- 1 - hazards
    
    # Calculate survival curves of all observed time points of each person
    predSurv <- aggregate(oneMinusPredHaz ~ obj, data = testDataLong, 
                          FUN=cumprod, na.action = NULL)
    
    # Calculate prediction error curve value in first interval
    pecObj <- predErrCurve(timepoints = 1, estSurvList = predSurv[[2]], 
                               testTime = testTime, testEvent = testEvent, 
                               trainTime = trainTime, trainEvent = trainEvent)
    
    # Estimate integrated prediction error curves
    if( is.null(tmax) ){
      tmax <- max(trainTime)
    }
    ipecVal <- unname(intPredErrCurve(predErrObj = pecObj, tmax = tmax))
    return(ipecVal)
}
