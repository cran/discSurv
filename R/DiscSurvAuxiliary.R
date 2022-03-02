#'Compute Subdistribution Weights 
#' 
#'@description Function to compute new subdistribution weights for a test data set based on the estimated 
#'censoring survival function from a learning data set 
#' 
#'@param dataShortTrain Learning data in short format ("class data.frame").
#'@param dataShortTest Test data in short format ("class data.frame"). 
#'@param timeColumn Character specifying the column name of the observed event times ("character vector"). 
#'It is required that the observed times are discrete ("integer vector").
#'@param eventColumns Character vector specifying the column names of the event indicators ("logical vector")(excluding censoring events). 
#'It is required that a 0-1 coding is used for all events. The algorithm treats row sums of zero of all event columns as censored.
#'@param eventFocus Column name of the event of interest, which corresponds to the type 1 event ("character vector").
#'@return Subdstribution weights for the test data in long format using the estimated 
#' censoring survival function from the learning data ("numeric vector"). The length of the 
#' vector is equal to the number of observations of the long test data. 
#'@author Moritz Berger \email{moritz.berger@@imbie.uni-bonn.de} \cr \url{https://www.imbie.uni-bonn.de/personen/dr-moritz-berger/}
#'@seealso \code{\link{dataLongSubDist}}, \code{\link{calibrationPlot}}
#'@references \insertRef{bergerSubdist}{discSurv}
#'@keywords discrete subdistribution hazards 
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
#' # Convert to long format
#' trainSet_long <- dataLongSubDist(dataShort = trainSet, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor4"), eventFocus = "censor1")
#' valSet_long <- dataLongSubDist(dataShort = valSet, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor4"), eventFocus = "censor1")
#' 
#' # Compute new weights of the validation data set 
#' valSet_long$subDistWeights <- weightsLtoT(trainSet, valSet, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor4"), eventFocus = "censor1")
#' 
#' # Estimate continuation ratio model with logit link
#' glmFit <- glm(formula = y ~ timeInt + age + logwage, data = trainSet_long, 
#' family = binomial(), weights = trainSet_long$subDistWeights)
#' 
#' # Calculate predicted discrete hazards 
#' predHazards <- predict(glmFit, newdata = valSet_long, type = "response")
#' 
#' # Calibration plot 
#' calibrationPlot(predHazards, testDataLong = valSet_long, weights = valSet_long$subDistWeights)
#' 
#' @export weightsLtoT
weightsLtoT <- function(dataShortTrain, dataShortTest, timeColumn, eventColumns, eventFocus){
  
  Ghat_l  <- estSurvCens(dataShort = dataShortTrain, timeColumn = timeColumn, eventColumns = eventColumns)
  
  data_t_timeColumn <- as.numeric(as.character(dataShortTest[, timeColumn]))
  tmax <- max(data_t_timeColumn)
  delta_i <- rowSums(dataShortTest[, eventColumns]) == 1
  epsilon_i <- sapply(1:dim(dataShortTest)[1], function(x) ifelse(delta_i[x], which(dataShortTest[x, eventColumns] == 1), 0))
  data_t_timeColumn <- ifelse(delta_i & dataShortTest[, eventFocus] == 0, tmax, data_t_timeColumn)
  obj <- rep(1:nrow(dataShortTest), each = tmax)
  timeInt <- c(sapply(1:length(data_t_timeColumn), function(k) 1:tmax))
  weights <- Ghat_l[timeInt]/Ghat_l[pmin(dataShortTest[obj, timeColumn], timeInt)] * (ifelse(timeInt <= dataShortTest[obj, timeColumn], 1, 0) + 
                                                                                 ifelse(dataShortTest[obj, timeColumn] <= (timeInt - 1) & delta_i[obj] * epsilon_i[obj] > 1, 1, 0))
  return(weights)
}
