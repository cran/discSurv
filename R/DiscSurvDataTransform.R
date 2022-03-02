#########################################################
# Transform continuous time variables to a discrete grid
#########################################################

# Description
# Discretizes continuous time variable into a specified grid of correct censored data

# Input
# dataShortTrain: Data frame with observed times, event indicator and covariates
# timeColumn: Character giving the column name of the observed times
# intervalLimits: Numeric vector of the right interval borders, e. g. if the intervals are
  # [0, a_1), [a_1, a_2), [a_2, a_{\max}), then intervalLimits = c(a_1, a_2, a_{\max})

# Output
# Original data frame with discretized response in the first column

#' Continuous to Discrete Transformation
#' 
#' Discretizes continuous time variable into a specified grid of censored data
#' for discrete survival analysis. It is a data preprocessing step, before the
#' data can be extendend in long format and further analysed with discrete
#' survival models.
#' 
#' 
#' @param dataShort Original data in short format ("class data.frame").
#' @param timeColumn Name of the column of discrete survival times ("character vector").
#' @param intervalLimits Right interval borders ("numeric vector"), e. g. if
#' the intervals are [0, a_1), [a_1, a_2), [a_2, a_max), then intervalLimits =
#' c(a_1, a_2, a_max)
#' @param equi Specifies if argument \emph{intervalLimits} should be interpreted as
#' number of equidistant intervals ("logical vector").
#' @param timeAsFactor Specifies if the computed discrete time intervals should be 
#' converted to a categorical variable ("logical vector"). Default is FALSE. 
#' In the default settings the discret time intervals are treated 
#' as quantitative ("numeric vector").
#' @return Gives the data set expanded with a first column "timeDisc". This
#' column includes the discrete time intervals ("class factor").
#' @note In discrete survival analysis the survival times have to be
#' categorized in time intervals. Therefore this function is required, if there
#' are observed continuous survival times.
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link{dataLong}}, \code{\link{dataLongTimeDep}},
#' \code{\link{dataLongCompRisks}}
#' @references
#' \insertRef{tutzModelDisc}{discSurv}
#' @keywords datagen
#' @examples
#' 
#' # Example copenhagen stroke study data
#' library(pec)
#' data(cost)
#' head(cost)
#' 
#' # Convert observed times to months
#' # Right borders of intervals [0, a_1), [a_1, a_2), ... , [a_{\\max-1}, a_{\\max})
#' IntBorders <- 1:ceiling(max(cost$time)/30)*30
#' 
#' # Select subsample
#' subCost <- cost [1:100, ]
#' CostMonths <- contToDisc(dataShort=subCost, timeColumn = "time", intervalLimits = IntBorders)
#' head(CostMonths)
#' 
#' # Select subsample giving number of equidistant intervals
#' CostMonths <- contToDisc(dataShort = subCost, timeColumn = "time", intervalLimits = 10, equi = TRUE)
#' head(CostMonths)
#' 
#' @export contToDisc
contToDisc <- function(dataShort, timeColumn, intervalLimits, equi=FALSE, timeAsFactor=FALSE) {
  
  # Input checks
  # Correct formats
  if(!is.data.frame(dataShort)) {stop("Argument *dataShort* is not in the correct format! Please specify as data.frame object.")}
  if(!(is.character(timeColumn))) {stop("Argument *timeColumn* is not in the correct format! Please specify as character.")}
  if(length(timeColumn) != 1) {stop("Argument *timeColumn* is not in the correct format! Please specify as scalar with length one.")}
  if(equi==FALSE) {
    if(!(is.vector(intervalLimits))) {stop("Argument *intervalLimits* is not in the correct format! Please specify as numeric vector.")}
    if(!(is.numeric(intervalLimits))) {stop("Argument *intervalLimits* is not in the correct format! Please specify as numeric vector.")}
  }
  if(equi == TRUE & length(intervalLimits) != 1) {stop("Argument *intervalLimits* is not in the correct format! Please specify a scalar integer value.")}
  if(equi == TRUE & any(floor(intervalLimits) != intervalLimits)) {stop("Argument *intervalLimits* is not in the correct format! Please specify an scalar integer value.")}
  
  # Can *timeColumn* be accessed in *dataShort*?
  if(any(class(tryCatch(dataShort [, timeColumn], error = function (e) e)) == "error"))
  {stop("*timeColumn* is not available in *dataShort*! Please specify the correct column of observed times.")}
 
  #####################
  # Main Code
  
  # Extract observed time
  obsTime <- dataShort [, timeColumn]
  
  if(equi == FALSE) {
    # Include zero in interval limits
    intervalLimits <- c(0, intervalLimits)
    
    # Check if breaks are not unique and add little bit of jitter in that case
    # Then the breaks are unique!
    Check <- tryCatch(cut(obsTime, intervalLimits, right = FALSE), error = function (e) "Bug!")
    if(identical(Check, "Bug!")) {
      intervalLimits [-1] <- jitter(intervalLimits [-1], factor = 0.01)
    }
  }
  
  # Divide the time to discrete intervals
  timeDisc <- as.numeric(cut(obsTime, intervalLimits, right = FALSE))
  if(timeAsFactor){
    timeDisc <- as.factor(timeDisc)
  }
 
  # Include discretized time in in original data set
  dataShort <- cbind(timeDisc = timeDisc, dataShort)
  return(dataShort)
}

##########################################################
# Data transformation of univeriate discrete survival data 
# without time varying covariates
##########################################################

# Description:
# Transform data in short format into long format for discrete survival analysis
# Data is assumed to include no time varying covariates, e. g. no follow up visits are allowed
# It is assumed that the covariates stay constant over time, in which no information is available

# Input
# dataShortTrain: Data frame with observed times, event indicator and covariates
# timeColumn: Character giving the column name of the observed times. 
  # It is required that the observed times are discrete (integer)
# eventColumn: Character giving the column name of the event indicator
 # It is required that this is a binary variable with 1=="event" and 0=="censored"

# Output: Original data.frame with three additional columns
# obj: Index of persons as integer vector
# timeInt: Index of time intervals as integer vector
# y: Response in long format as binary vector. 1=="event happens in period timeInt" and zero otherwise



#' Data Long Transformation
#' 
#' Transform data from short format into long format for discrete survival
#' analysis and right censoring. Data is assumed to include no time varying
#' covariates, e. g. no follow up visits are allowed. It is assumed that the
#' covariates stay constant over time, in which no information is available.
#' 
#' If the data has continuous survival times, the response may be transformed
#' to discrete intervals using function \code{\link{contToDisc}}. If the data
#' set has time varying covariates the function \code{\link{dataLongTimeDep}}
#' should be used instead. In the case of competing risks and no time varying
#' covariates see function \code{\link{dataLongCompRisks}}.
#' 
#' @param dataShort Original data in short format ("class data.frame").
#' @param timeColumn Character giving the column name of the observed times. It
#' is required that the observed times are discrete ("integer vector").
#' @param eventColumn Column name of the event indicator ("character vector").
#' It is required that this is a binary variable with 1=="event" and
#' 0=="censored".
#' @param timeAsFactor Should the time intervals be coded as factor ("logical vector")? 
#' Default is FALSE. In the default settings the column is treated as quantitative variable ("numeric vector").
#' @param remLastInt Should the last theoretical interval be removed in long
#' format ("logical vector")? Default setting (FALSE) is no deletion. This is only important, if the short format
#' data includes the last theoretic interval [a_q, Inf). There are only events
#' in the last theoretic interval, so the discrete hazard is always one and these
#' observations have to be excluded for estimation.
#' @param aggTimeFormat Instead of the usual long format, should every
#' observation have all time intervals ("logical vector")? Default is standard
#' long format (FALSE). In the case of nonlinear risk score models, the time effect has
#' to be integrated out before these can be applied to the C-index.
#' @param lastTheoInt Gives the number of the last theoretic interval ("integer vector"). 
#' Only used, if argument \emph{aggTimeFormat} is set to TRUE.
#' @return Original data.frame with three additional columns: \itemize{ \item
#' {obj} Index of persons as integer vector \item {timeInt} Index of time
#' intervals (factor) \item {y} Response in long format as binary vector.
#' 1=="event happens in period timeInt" and zero otherwise. 
#' If argument \emph{responseAsFactor} is set to TRUE, then responses will be coded as factor 
#' in one column.}
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' 
#' Matthias Schmid \email{matthias.schmid@@imbie.uni-bonn.de}
#' @seealso \code{\link{contToDisc}}, \code{\link{dataLongTimeDep}},
#' \code{\link{dataLongCompRisks}}
#' @references
#' \insertRef{tutzModelDisc}{discSurv} \cr\cr
#' \insertRef{fahrmeirDiscSurv}{discSurv} \cr\cr
#' \insertRef{thompsonTreatment}{discSurv}
#' @keywords datagen
#' @examples
#' 
#' # Example unemployment data
#' library(Ecdat)
#' data(UnempDur)
#' 
#' # Select subsample
#' subUnempDur <- UnempDur [1:100, ]
#' head(subUnempDur)
#' 
#' # Convert to long format
#' UnempLong <- dataLong (dataShort = subUnempDur, timeColumn = "spell", eventColumn = "censor1")
#' head(UnempLong, 20)
#' 
#' # Is there exactly one observed event of y for each person?
#' splitUnempLong <- split(UnempLong, UnempLong$obj)
#' all(sapply(splitUnempLong, function (x) sum(x$y))==subUnempDur$censor1) # TRUE
#' 
#' # Second example: Acute Myelogenous Leukemia survival data
#' library(survival)
#' head(leukemia)
#' leukLong <- dataLong(dataShort = leukemia, timeColumn = "time", 
#' eventColumn = "status", timeAsFactor=TRUE)
#' head(leukLong, 30)
#' 
#' # Estimate discrete survival model
#' estGlm <- glm(formula = y ~ timeInt + x, data=leukLong, family = binomial())
#' summary(estGlm)
#' 
#' # Estimate survival curves for non-maintained chemotherapy
#' newDataNonMaintained <- data.frame(timeInt = factor(1:161), x = rep("Nonmaintained"))
#' predHazNonMain <- predict(estGlm, newdata = newDataNonMaintained, type = "response")
#' predSurvNonMain <- cumprod(1-predHazNonMain)
#' 
#' # Estimate survival curves for maintained chemotherapy
#' newDataMaintained <- data.frame(timeInt = factor(1:161), x = rep("Maintained"))
#' predHazMain <- predict(estGlm, newdata = newDataMaintained, type = "response")
#' predSurvMain <- cumprod(1-predHazMain)
#' 
#' # Compare survival curves
#' plot(x = 1:50, y = predSurvMain [1:50], xlab = "Time", ylab = "S(t)", las = 1, 
#' type = "l", main = "Effect of maintained chemotherapy on survival of leukemia patients")
#' lines(x = 1:161, y = predSurvNonMain, col = "red")
#' legend("topright", legend = c("Maintained chemotherapy", "Non-maintained chemotherapy"), 
#' col = c("black", "red"), lty = rep(1, 2))
#' # The maintained therapy has clearly a positive effect on survival over the time range
#' 
#' ##############################################
#'# Simulation
#'# Single event in case of right-censoring
#'
#'# Simulate multivariate normal distribution
#'library(discSurv)
#'library(mvnfast)
#'set.seed(-1980)
#'X <- mvnfast::rmvn(n = 1000, mu = rep(0, 10), sigma = diag(10))
#'
#'
#'# Specification of discrete hazards with 11 theoretical intervals
#'betaCoef <- seq(-1, 1, length.out = 11)[-6]
#'timeInt <- seq(-1, 1, length.out = 10)
#'linPred <- c(X %*% betaCoef)
#'hazTimeX <- cbind(sapply(1:length(timeInt), 
#'                         function(x) exp(linPred+timeInt[x]) / (1+exp(linPred+timeInt[x])) ), 1)
#'
#'
#'# Simulate discrete survival and censoring times in 10 observed intervals
#'discT <- rep(NA, dim(hazTimeX)[1])
#'discC <- rep(NA, dim(hazTimeX)[1])
#'for( i in 1:dim(hazTimeX)[1] ){
#'  
#'  discT[i] <- sample(1:11, size = 1, prob = estMargProb(haz=hazTimeX[i, ]))
#'  discC[i] <- sample(1:11, size = 1, prob = c(rep(1/11, 11)))
#'}
#'
#'
#'# Calculate observed times, event indicator and specify short data format
#'eventInd <- discT <= discC
#'obsT <- ifelse(eventInd, discT, discC)
#'eventInd[obsT == 11] <- 0
#'obsT[obsT == 11] <- 10
#'simDatShort <- data.frame(obsT = obsT, event = as.numeric(eventInd), X)
#'
#'
#'# Convert data to discrete data long format
#'simDatLong <- dataLong(dataShort = simDatShort, timeColumn = "obsT", eventColumn = "event",
#'timeAsFactor=TRUE)
#'
#'
#'# Estimate discrete-time continuation ratio model
#'formSpec <- as.formula(paste("y ~ timeInt + ", 
#'                             paste(paste("X", 1:10, sep=""), collapse = " + "), sep = ""))
#'modelFit <- glm(formula = formSpec, data = simDatLong, family = binomial(link = "logit"))
#'summary(modelFit)
#'
#'
#'# Compare estimated to true coefficients
#'coefModel <- coef(modelFit)
#'MSE_covariates <- mean((coefModel[11:20]-timeInt)^2)
#'MSE_covariates
#' # -> Estimated coefficients are near true coefficients
#' 
#' @export dataLong
dataLong <- function(dataShort, timeColumn, eventColumn, timeAsFactor=FALSE,
                      remLastInt=FALSE, aggTimeFormat=FALSE, lastTheoInt=NULL) {

  # Input checks
  if(!is.data.frame(dataShort)) {stop("Argument *dataShort* is not in the correct format! Please specify as data.frame object.")}
  if(!(is.character(timeColumn))) {stop("Argument *timeColumn* is not in the correct format! Please specify as character.")}
  if(length(timeColumn)!=1) {stop("Argument *timeColumn* is not in the correct format! Please specify as scalar with length one.")}
  if(!(is.character(eventColumn))) {stop("Argument *eventColumn* is not in the correct format! Please specify as character.")}
  if(length(eventColumn)!=1) {stop("Argument *eventColumn* is not in the correct format! Please specify as scalar with length one.")}

  # Can *timeColumn* be accessed in *dataShort*?
  if(any(class(tryCatch(dataShort [, timeColumn], error=function (e) e)) == "error")) {
    stop("*timeColumn* is not available in *dataShort*! Please specify the correct column of observed times.")
  }
  
  # Can *eventColumn* be accessed in *dataShort*?
  if(any(class(tryCatch(dataShort [, eventColumn], error = function (e) e)) == "error")) {
    stop("*eventColumn* is not available in *dataShort*! Please specify the correct column of the event indicator.")
  }
  
  # Check if observed times are only integer values
  if(!all(dataShort [, timeColumn] == floor(as.numeric(as.character(dataShort [, timeColumn]))))) {
    stop("*timeColumn* has not only integer values! Please convert the observed time in discrete intervals.")
  }
  
  # Check if event indicator has only zero or one values
  if(!all(as.numeric(as.character(dataShort [, eventColumn])) == 0 | as.numeric(as.character(dataShort [, eventColumn])) == 1)) {
    stop("*eventColumn* is not a binary vector! Please check, that events equals 1 and 0 otherwise.")
  }
  
  ##################
  # Main Code
  
  # Extract column index of survival and censoring times
  c1 <- which(eval(timeColumn) == colnames(dataShort))
  c2 <- which(eval(eventColumn) == colnames(dataShort))
  
  # Construct indices of persons
  if(aggTimeFormat){
    obj <- rep(1:nrow(dataShort), each = lastTheoInt)
  } else{
    obj <- rep(1:nrow(dataShort), as.vector(dataShort[, c1]))
  }

  # Long format of covariates
  dataShortLong <- dataShort[obj, ]
  
  # Calculate discrete time interval
  if(aggTimeFormat){
    if(timeAsFactor) {
      timeInt <- factor( rep(1:lastTheoInt, nrow(dataShort)) )
    }
    else{
      timeInt <- rep(1:lastTheoInt, nrow(dataShort))
    }
    # Calculate response
    y <- c(unlist(apply(dataShort, 1, FUN = 
        function(k) c(rep(0, as.numeric(k[c1])-1), as.numeric(k[c2]), 
                      rep(0, lastTheoInt - as.numeric(k[c1]) ) ) )))
    
  } else{
    if(timeAsFactor) {
      timeInt <- factor(unlist(apply(dataShort, 1, FUN = function(k) 1:k[c1])))
    }
    else{
      timeInt <- unlist(apply(dataShort, 1, FUN = function(k) 1:k[c1]))
    }
    # Calculate response
    y <- c(unlist(apply(dataShort, 1, FUN = 
         function(k) c(rep(0, as.numeric(k[c1])-1), as.numeric(k[c2])) )))
  }
  
  # Aggregate results in one data.frame
  dataShortLong <- cbind(obj, timeInt, y, dataShortLong)
  
  if(remLastInt) {
    # Remove cases with observed values in last interval, 
    # because the discrete hazard is always 1
    # t == t_max & event == 1 -> t == t_max-1 & event == 0
    remInd <- which(dataShortLong$y == 1 & 
                      dataShortLong$timeInt == max(as.numeric(as.character(
                        dataShortLong$timeInt))) )
    if(length(remInd) != 0) {
      dataShortLong <- dataShortLong[-remInd, ]
    }
  }
  
  return(dataShortLong)
} 

#########################################################
# Datalong for time dependent covariates
#########################################################

# Description:
# Transform data in short format into long format for discrete survival analysis
# Data is assumed to include time varying covariates, e. g. follow up visits are allowed
# It is assumed that the covariates stay constant between measurements, where no information is available

# Input
# dataShortTrain: Data frame with observed times, event indicator and covariates
# timeColumn: Character giving the column name of the observed times. 
# It is required that the observed times are discrete.
# eventColumn: Character giving the column name of the event indicator
# It is required that this is a binary variable with 1=="event" and 0=="censored"
# idColumn: Gives the identification number of the persons

# Output: Original data.frame with three additional columns
# obj: Index of persons as integer vector
# timeInt: Index of time intervals as integer vector
# y: Response in long format as binary vector. 1 == "event happens in period timeInt" and 0 otherwise



#' Data Long Time Dependent Covariates
#' 
#' Transforms short data format to long format for discrete survival modelling
#' of single event analysis with right censoring. Covariates may vary over
#' time.
#' 
#' There may be some intervals, where no additional information on the
#' covariates is observed (e. g. observed values in interval one and three but
#' two is missing). In this case it is assumed, that the values from the last
#' observation stay constant over time until a new measurement was done.
#' 
#' @param dataSemiLong Original data in semi-long format ("class data.frame").
#' @param timeColumn Character giving the column name of the observed times ("character vector"). 
#' It is required that the observed times are discrete ("integer vector").
#' @param eventColumn Column name of the event indicator ("character vector").
#' It is required that this is a binary variable with 1=="event" and
#' 0=="censored".
#' @param idColumn Name of column of identification number of persons ("character vector").
#' @param timeAsFactor Should the time intervals be coded as factor ("logical vector")? Default is
#' FALSE. In case of default settings the discrete time intervals are treated as quantitative ("numeric vector"). 
#' @return Original data in long format with three additional columns: \itemize{ \item
#' {obj} Index of persons as integer vector \item {timeInt} Index of time
#' intervals (factor) \item {y} Response in long format as binary vector.
#' 1=="event happens in period timeInt" and zero otherwise }
#' @details In contrast to continuous survival (see e. g. \code{\link[survival]{Surv}}) 
#' the start and stop time notation is not used here. In discrete time survival analysis the only relevant
#' information is to use the stop time. Start time does not matter, because all discrete intervals need to be  
#' included in the long data set format to ensure consistent estimation. It is assumed that the supplied 
#' data set "dataSemiLong" contains all repeated measurements of each cluster in semi-long format (e. g. persons). 
#' For further information see example \emph{Start-stop notation}.
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link{contToDisc}}, \code{\link{dataLong}},
#' \code{\link{dataLongCompRisks}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv} \cr\cr
#' \insertRef{fahrmeirDiscSurv}{discSurv} \cr\cr
#' \insertRef{thompsonTreatment}{discSurv}
#' @keywords datagen
#' @examples
#' 
#' # Example Primary Biliary Cirrhosis data
#' library(survival)
#' dataSet1 <- pbcseq
#' 
#' # Only event death is of interest
#' dataSet1$status [dataSet1$status == 1] <- 0
#' dataSet1$status [dataSet1$status == 2] <- 1
#' table(dataSet1$status)
#' 
#' # Convert to months
#' dataSet1$day <- ceiling(dataSet1$day/30) + 1
#' names(dataSet1) [7] <- "month"
#' 
#' # Convert to long format for time varying effects
#' pbcseqLong <- dataLongTimeDep (dataSemiLong = dataSet1, timeColumn = "month", 
#' eventColumn = "status", idColumn = "id")
#' pbcseqLong [pbcseqLong$obj == 1, ]
#' 
#' #####################
#' # Start-stop notation
#' 
#' library(survival)
#' ?survival::heart
#' 
#' # Assume that time was measured on a discrete scale.
#' # Discrete interval lengths are assumed to vary.
#' intervalLimits <- quantile(heart$stop, probs = seq(0.1, 1, by=0.1))
#' intervalLimits[length(intervalLimits)] <- intervalLimits[length(intervalLimits)] + 1
#' heart_disc <- contToDisc(dataShort = heart, timeColumn = "stop", 
#' intervalLimits = intervalLimits, equi = FALSE)
#' table(heart_disc$timeDisc)
#' 
#' # Conversion to long format
#' heart_disc_long <- dataLongTimeDep(dataSemiLong = heart_disc, timeColumn = "timeDisc", 
#' eventColumn = "event", idColumn = "id")
#' head(heart_disc_long)
#' 
#' @export dataLongTimeDep
dataLongTimeDep <- function(dataSemiLong, timeColumn, eventColumn, idColumn, timeAsFactor=FALSE) {

  # Input checks
  if(!is.data.frame(dataSemiLong)) {stop("Argument *dataSemiLong* is not in the correct format! Please specify as data.frame object.")}
  if(!(is.character(timeColumn))) {stop("Argument *timeColumn* is not in the correct format! Please specify as character.")}
  if(length(timeColumn)!=1) {stop("Argument *timeColumn* is not in the correct format! Please specify as scalar with length one.")}
  if(!(is.character(eventColumn))) {stop("Argument *eventColumn* is not in the correct format! Please specify as character.")}
  if(length(eventColumn)!=1) {stop("Argument *eventColumn* is not in the correct format! Please specify as scalar with length one.")}
  
  # Can *timeColumn* be accessed in *dataSemiLong*?
  if(any(class(tryCatch(dataSemiLong [, timeColumn], error=function (e) e)) == "error")) {
    stop("*timeColumn* is not available in *dataSemiLong*! Please specify the correct column of observed times.")
  }
  # Can *eventColumn* be accessed in *dataSemiLong*?
  if(any(class(tryCatch(dataSemiLong [, eventColumn], error=function (e) e)) == "error")) {
    stop("*eventColumn* is not available in *dataSemiLong*! Please specify the correct column of the event indicator.")
  }
  # Can *idColumn* be accessed in *dataSemiLong*?
  if(any(class(tryCatch(dataSemiLong [, idColumn], error = function (e) e)) == "error")) {
    stop("*idColumn* is not available in *dataSemiLong*! Please specify the correct column of the identification number.")
  }
  
  # Check if observed times are only integer values
  if(!all(dataSemiLong [, timeColumn] == floor(as.numeric(as.character(dataSemiLong [, timeColumn]))))) {
    stop("*timeColumn* has not only integer values! Please convert the observed time in discrete intervals.")
  }
  # Check if event indicator has only zero or one values
  if(!all(as.numeric(as.character(dataSemiLong [, eventColumn])) == 0 | as.numeric(as.character(dataSemiLong [, eventColumn])) == 1)) {
    stop("*eventColumn* is not a binary vector! Please check, that events equals 1 and 0 otherwise.")
  }
  
  ###############
  # Main code

  # Rearrange original data, that ID is in increasing order
  dataSemiLong <- dataSemiLong[order(dataSemiLong [, idColumn]), ]
  
  # Split data by persons
  splitDataSet <- split(x = dataSemiLong, f = dataSemiLong [, idColumn])
  lengthSplitDataSet <- 1:length(splitDataSet)
  
  # Get count of replicates in each split
  Counts <- lapply(splitDataSet, function (x) c(diff(x [, timeColumn]), 1))
  SumCounts <- as.numeric(sapply(Counts, function (x) sum(x)))
  
  # Replicate indices in each split
  Indizes <- lapply(lengthSplitDataSet, function (x) rep(x = 1:dim(splitDataSet [[x]])[1], times = Counts [[x]]))
  
  # Duplicate rows for each split and combine the results
  dataList <- lapply(lengthSplitDataSet, function (x) splitDataSet [[x]] [Indizes [[x]], ])
  dataList <- do.call(rbind, dataList)
  
  # Create ID variable in long format for each split and combine the results
  dataListID <- lapply(lengthSplitDataSet, function (x) rep(x, SumCounts [x]))
  dataListID <- do.call(c, dataListID)
  
  # Create time interval variable in long format for each split and combine results
  dataListTimeInt <- lapply(lengthSplitDataSet, function (x) 1:SumCounts [x])
  if(timeAsFactor) {
    dataListTimeInt <- factor(do.call(c, dataListTimeInt))
  }
  else{
    dataListTimeInt <- do.call(c, dataListTimeInt)
  }
  
  # Create time response variable in long format for each split and combine results
  dataListResponse <- lapply(lengthSplitDataSet, function (x) c(rep(0, SumCounts [x]-1), as.numeric(as.character(tail(splitDataSet [[x]] [, eventColumn], 1)))))
  dataListResponse <- do.call(c, dataListResponse)
  
  # Combine overall results and output
  dataComplete <- cbind(obj = dataListID, timeInt = dataListTimeInt, y = dataListResponse, dataList)
  return(dataComplete)
}

##############################################################
# Function for dataLong in the case of competing risks models
##############################################################

# Description
# Constructs from a short data format a long format for discrete survival modelling
# Assumptions: 
# 1. Censoring process is independent of survival process
# 2. non informative censoring
# 3. correct censoring
# 4. covariates dont vary over time

# Input
# dataShortTrain: Complete data set in short format (data.frame)
# timeColumn: Character scalar of column names of numeric time variable
# eventColumns: Character vector of column names of event indicators (without indicator for censoring)
# Baseline (all event indicators are zero) is assumend to be interpreted as censored observation

# Output
# Data set in long format (data.frame)
# obj: Gives identification number of objects (row index in short format) (integer)
# timeInt: Gives number of discrete time interval (factor)
# e0: No event (observation censored in specific interval)
# e1: Indicator of first event, 1 if event takes place and 0 otherwise
# ...
# ei: Indicator of i-th event, 1 if event takes place and 0 otherwise
# ...
# ek: Indicator of last event, 1 if event takes place and 0 otherwise
# dataShortTrain: Original data with replicated rows



#' Data Long Competing Risks Transformation
#' 
#' Transforms short data format to long format for discrete survival modelling
#' in the case of competing risks with right censoring. It is assumed that the
#' covariates are not time varying.
#' 
#' It is assumed, that only one event happens at a specific time point
#' (competing risks). Either the observation is censored or one of the possible
#' events takes place.
#' 
#' @param dataShort Original data in short format ("class data.frame").
#' @param timeColumn Character giving the column name of the observed times ("character vector"). It
#' is required that the observed times are discrete ("integer vector").
#' @param eventColumns Character vector giving the column names of the event
#' indicators (excluding censoring column)("character vector"). It is required that all events are
#' binary encoded. If the sum of all event indicators is zero, then this is
#' interpreted as a censored observation. Alternatively a column name of a
#' factor representing competing events can be given. In this case the argument
#' \emph{eventColumnsAsFactor} has to be set TRUE and the first level is assumed to
#' represent censoring.
#' @param eventColumnsAsFactor Should the argument \emph{eventColumns} be interpreted
#' as column name of a factor variable ("logical vector")? Default is FALSE.
#' @param timeAsFactor Should the time intervals be coded as factor ("logical vector")? 
#' Default is FALSE. In the default settings the discrete time variable are treated as quantitative.
#' @param aggTimeFormat Instead of the usual long format, should every
#' observation have all time intervals ("logical vector")? Default is standard
#' long format. In the case of nonlinear risk score models, the time effect has
#' to be integrated out before these can be applied to the C-index.
#' @param lastTheoInt Gives the number of the last theoretic interval ("integer
#' vector"). Only used, if \emph{aggTimeFormat} is set to TRUE.
#' @param responseAsFactor Should the response columns be given as factor ("logical vector")? 
#' Default is FALSE.
#' @return Original data set in long format with additional columns \itemize{
#' \item {obj} Gives identification number of objects (row index in short
#' format) (integer) \item {timeInt} Gives number of discrete time intervals
#' (factor) \item {responses} Columns with dimension count of events + 1
#' (censoring) \itemize{ \item {e0} No event (observation censored in specific
#' interval) \item {e1} Indicator of first event, 1 if event takes place and 0
#' otherwise \item ... ...  \item {ek} Indicator of last k-th event, 1 if event
#' takes place and zero otherwise}
#' If argument responseAsFactor=TRUE, then responses will be coded as factor in one column.
#' }
#' @details In contrast to continuous survival (see e. g. \code{\link[survival]{Surv}}) 
#' the start and stop time notation is not used here. In discrete time survival analysis the only relevant
#' information is to use the stop time. Start time does not matter, because all discrete intervals need to be  
#' included in the long data set format to ensure consistent estimation. It is assumed that the supplied 
#' data set \emph{dataShort} contains all repeated measurements of each cluster (e. g. persons). 
#' For further information see example \emph{Start-stop notation}.
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link{contToDisc}}, \code{\link{dataLongTimeDep}},
#' \code{\link{dataLongCompRisksTimeDep}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv} \cr\cr
#' \insertRef{steeleMultistate}{discSurv} \cr\cr
#' \insertRef{wijiUnemployment}{discSurv}
#' @keywords datagen
#' @examples
#' 
#' # Example with unemployment data
#' library(Ecdat)
#' data(UnempDur)
#' 
#' # Select subsample
#' SubUnempDur <- UnempDur [1:100, ]
#' 
#' # Convert competing risk data to long format
#' SubUnempDurLong <- dataLongCompRisks (dataShort = SubUnempDur, timeColumn = "spell", 
#' eventColumns = c("censor1", "censor2", "censor3", "censor4"))
#' head(SubUnempDurLong, 20)
#' 
#' # Fit multinomial logit model with VGAM package
#' # with one coefficient per response
#' library(VGAM)
#' multLogitVGM <- vgam(cbind(e0, e1, e2, e3, e4) ~ timeInt + ui + age + logwage,
#'                     family = multinomial(refLevel = 1), 
#'                     data = SubUnempDurLong)
#' coef(multLogitVGM)
#' 
#' # Alternative: Use nnet
#' # Convert response to factor
#' rawResponseMat <- SubUnempDurLong[, c("e0", "e1", "e2", "e3", "e4")]
#' NewFactor <- factor(unname(apply(rawResponseMat, 1, function(x) which(x == 1))), 
#'                     labels = colnames(rawResponseMat))
#' 
#' # Include recoded response in data
#' SubUnempDurLong <- cbind(SubUnempDurLong, NewResp = NewFactor)
#' 
#' # Construct formula of mlogit model
#' mlogitFormula <- formula(NewResp ~ timeInt + ui + age + logwage)
#' 
#' # Fit multinomial logit model
#' # with one coefficient per response
#' library(nnet)
#' multLogitNNET <- multinom(formula = mlogitFormula, data = SubUnempDurLong)
#' coef(multLogitNNET)
#' 
#' ###########################################################
#'# Simulation
#'# Cause specific competing risks in case of right-censoring
#'# Discrete subdistribution hazards model
#'
#'# Simulate covariates as multivariate normal distribution
#'library(mvnfast)
#'set.seed(1980)
#'X <- mvnfast::rmvn(n = 1000, mu = rep(0, 4), sigma = diag(4))
#'
#' # Specification of two discrete cause specific hazards with four intervals
#' # Event 1
#'theoInterval <- 4
#'betaCoef_event1 <- seq(-1, 1, length.out = 5)[-3]
#'timeInt_event1 <- seq(0.1, -0.1, length.out = theoInterval-1)
#'linPred_event1 <- c(X %*% betaCoef_event1)
#'# Event 2
#'betaCoef_event2 <- seq(-0.5, 0.5, length.out = 5)[-3]
#'timeInt_event2 <- seq(-0.1, 0.1, length.out = theoInterval-1)
#'linPred_event2 <- c(X %*% betaCoef_event2)
#'# Discrete cause specific hazards in last theoretical interval
#'theoHaz_event1 <- 0.5
#'theoHaz_event2 <- 0.5
#'
# Derive discrete all cause hazard
#'haz_event1_X <- cbind(sapply(1:length(timeInt_event1), 
#'                function(x) exp(linPred_event1 + timeInt_event1[x]) / 
#'                (1 + exp(linPred_event1 + timeInt_event1[x]) + 
#'                exp(linPred_event2 + timeInt_event2[x])) ), theoHaz_event1)
#'
#'haz_event2_X <- cbind(sapply(1:length(timeInt_event2), 
#'                function(x) exp(linPred_event2 + timeInt_event2[x]) / 
#'                (1 + exp(linPred_event1 + timeInt_event1[x]) + 
#'                exp(linPred_event2 + timeInt_event2[x]) ) ), theoHaz_event2)
#'allCauseHaz_X <- haz_event1_X + haz_event2_X
#'
#'
# Calculate all cause probability P(T = t | X)
#'pT_X <- t(sapply(1:dim(allCauseHaz_X)[1], function(i) estMargProb(allCauseHaz_X[i, ]) ))
#'
#'
# Calculate event probability given time interval P(R = r | T = t, X)
#'pR_T_X_event1 <- haz_event1_X / (haz_event1_X + haz_event2_X)
#'
#'
# Simulate discrete survival times
#'survT <- sapply(1:dim(pT_X)[1], function(i) sample(x = 1:(length(timeInt_event1) + 1), 
#'                                                   size = 1, prob = pT_X[i, ]) )
#'censT <- sample(x = 1:(length(timeInt_event1)+1), size = dim(pT_X)[1], 
#'                prob = rep(1/(length(timeInt_event1) + 1), (length(timeInt_event1) + 1)), 
#'                replace = TRUE)
#'
#'
# Calculate observed times
#'obsT <- ifelse(survT <= censT, survT, censT)
#'obsEvent <- rep(0, length(obsT))
#'obsEvent <- sapply(1:length(obsT), 
#'                   function(i) if(survT[i] <= censT[i]){
#'                     return(sample(x = c(1, 2), size=1, 
#'                     prob = c(pR_T_X_event1[i, obsT[i]  ], 
#'                     1 - pR_T_X_event1[i, obsT[i]  ]) ))
#'                   } else{
#'                     
#'                     return(0)
#'                   }
#')
#'
#'
#'# Recode last interval to censored
#'lastInterval <- obsT == theoInterval
#'obsT[lastInterval] <- theoInterval - 1
#'obsEvent[lastInterval] <- 0
#'obsT <- factor(obsT)
#'obsEvent <- factor(obsEvent)
#'
# Conversion to long data format
#'datShort <- data.frame(event = factor(obsEvent), time = obsT, X)
#'datLong <- dataLongCompRisks(dataShort = datShort, timeColumn = "time", 
#'                             eventColumns = "event", responseAsFactor = TRUE, 
#'                             eventColumnsAsFactor = TRUE, timeAsFactor = TRUE)
#'
#'
#'# Estimate discrete cause specific hazard model
#'library(VGAM)
#'estModel <- vglm(formula=responses ~ timeInt + X1 + X2 + X3 + X4, data=datLong, 
#'                 family = multinomial(refLevel = 1))
#'
#'
#' # Mean squared errors per event
#' coefModels <- coef(estModel)
#' mean((coefModels[seq(7, length(coefModels), 2)] - betaCoef_event1)^2) # Event 1
#' mean((coefModels[seq(8, length(coefModels), 2)] - betaCoef_event2)^2) # Event 2
#' # -> Estimated coefficients are near true coefficients for each event type
#' 
#' @export dataLongCompRisks
dataLongCompRisks <- function (dataShort, timeColumn, eventColumns, eventColumnsAsFactor = FALSE, 
                                timeAsFactor = FALSE, aggTimeFormat = FALSE, lastTheoInt = NULL, 
                                responseAsFactor = FALSE) {
  
  # Input checks
  if(!is.data.frame(dataShort)) {stop("Argument *dataShort* is not in the correct format! Please specify as data.frame object.")}
  if(!(is.character(timeColumn))) {stop("Argument *timeColumn* is not in the correct format! Please specify as character.")}
  if(length(timeColumn)!=1) {stop("Argument *timeColumn* is not in the correct format! Please specify as scalar with length one.")}
  if(!all(is.character(eventColumns))) {stop("Argument *eventColumns* is not in the correct format! Please specify as character.")}
  
  # Can *timeColumn* be accessed in *dataShort*?
  if(any(class(tryCatch(dataShort [, timeColumn], error=function (e) e)) == "error")) {
    stop("*timeColumn* is not available in *dataShort*! Please specify the correct column of observed times.")
  }
  # Can *eventColumns* be accessed in *dataShort*?
  if(any(class(tryCatch(dataShort [, eventColumns], error=function (e) e)) == "error")) {
    stop("*eventColumns* is not available in *dataShort*! Please specify the correct column of the event indicator.")
  }
  
  # Check if observed times are only integer values
  if(!all(dataShort [, timeColumn]==floor(as.numeric(as.character(dataShort [, timeColumn]))))) {
    stop("*timeColumn* has not only integer values! Please convert the observed time in discrete intervals.")
  }
  
  #############
  # Main Code
  
  # Alternative input for competing risks coding 
  # as factor instead of multiple binary columns
  if(eventColumnsAsFactor) {
    respFact <- dataShort[, eventColumns]
    respMat <- model.matrix(~., data.frame(respFact))[, -1]
    eventColumns <- levels(respFact)[-1]
    dimnames(respMat) [[2]] <- eventColumns
    dataShort <- cbind(respMat, dataShort)
  }
  
  # Check if event indicators have only zero or one values
  checkVec <- vector("logical", length(eventColumns))
  for(i in 1:length(eventColumns)) {
    checkVec [i] <- all(as.numeric(as.character(dataShort [, eventColumns [i] ])) == 0 | as.numeric(as.character(dataShort [, eventColumns [i] ])) == 1)
  }
  if(!all(checkVec)) {
    stop("*eventColumns* is not a binary vector! Please check, that events equals 1 and 0 otherwise.")
  }
  
  # Index of columns of timeColumn and event indicators
  indextimeColumn <- which(eval(timeColumn) == colnames(dataShort))
  indizeseventColumns <- sapply(1:length(eventColumns), function (x) which(eval(eventColumns [x]) == colnames(dataShort)))
  dataSet_timeColumn <- as.numeric(as.character(dataShort[, indextimeColumn]))
  
  # Construct object counter, covariates in long format and discrete time intervals
  if(aggTimeFormat){
    
    obj <- rep(1:nrow(dataShort), each = lastTheoInt)
    dataSetLong <- dataShort[obj,]
    timeInt <- rep(1:lastTheoInt, nrow(dataShort))
    if(timeAsFactor){timeInt <- factor(timeInt)}
    
  } else{
    
    obj <- rep(1:nrow(dataShort), dataSet_timeColumn)
    dataSetLong <- dataShort[obj,]
    timeInt <- unlist(apply(dataShort, 1, FUN = function(k) 1:k[indextimeColumn]))
    if(timeAsFactor){timeInt <- factor(timeInt)}
    
  }
  
  # Construct censoring variable in short format
  dataSet_eventColumns <- dataShort [, indizeseventColumns]
  eventColumnsorShort <- 1 - rowSums(dataSet_eventColumns)
  responsesShort <- as.matrix(cbind(eventColumnsorShort, dataSet_eventColumns))
  dimnames(responsesShort) [[2]] <- 1:length(dimnames(responsesShort) [[2]])
  
  # Construct responses and censoring variable
  if(aggTimeFormat){
    
    NoeventColumns <- length(indizeseventColumns)
    responses <- matrix(0, nrow = 0, ncol = NoeventColumns + 1)
    for(i in 1:dim(dataShort) [1]) {
      row_rep <- as.numeric(as.character(dataShort [i, indextimeColumn])) - 1
      mat_temp <- matrix(rep(c(1, rep(0, NoeventColumns)), row_rep), nrow = row_rep, ncol = NoeventColumns + 1, byrow = TRUE)
      mat_temp <- rbind(mat_temp, responsesShort [i, ])
      mat_temp_after <- matrix(rep(c(1, rep(0, NoeventColumns)), lastTheoInt-row_rep-1), 
                               nrow = lastTheoInt-row_rep-1, ncol = NoeventColumns + 1, byrow = TRUE)
      responses <- rbind(responses, mat_temp, mat_temp_after)
    }
    dimnames(responses) [[2]] <- paste("e", 0:NoeventColumns, sep = "")
    
  } else{
    
    NoeventColumns <- length(indizeseventColumns)
    responses <- matrix(0, nrow = 0, ncol = NoeventColumns + 1)
    for(i in 1:dim(dataShort) [1]) {
      row_rep <- as.numeric(as.character(dataShort [i, indextimeColumn])) - 1
      mat_temp <- matrix(rep(c(1, rep(0, NoeventColumns)), row_rep), nrow = row_rep, ncol = NoeventColumns + 1, byrow = TRUE)
      mat_temp <- rbind(mat_temp, responsesShort [i, ])
      responses <- rbind(responses, mat_temp)
    }
    dimnames(responses) [[2]] <- paste("e", 0:NoeventColumns, sep = "")
    
  }
  
  # Remove unnecessary columns
  if(eventColumnsAsFactor) {
    dataSetLong <- dataSetLong[, -which(names(dataSetLong) %in% eventColumns)]
  }
  
  # Transform response columns to factor
  if(responseAsFactor) {
    colName <- colnames(responses)
    colNameInd <- sapply(1:dim(responses)[1], 
                         function(x) which(responses[x, ] == 1) )
    responses <- factor(colName[colNameInd], levels = colName)
  }
  
  # Combine results and output
  dataSetLong <- cbind(obj, timeInt, responses, dataSetLong)
  return(dataSetLong)  
}


##################################################################################
# Function for dataLong in the case of competing risks models 
# with time dependent covariates
##################################################################################

#' Data Long Competing Risks Time Dependent Covariates Transformation
#' 
#' Transforms short data format to long format for discrete survival modelling
#' in the case of competing risks with right censoring. Covariates may vary
#' over time.
#' 
#' There may be some intervals, where no additional information on the
#' covariates is observed (e. g. observed values in interval one and three but
#' two is missing). In this case it is assumed, that the values from the last
#' observation stay constant over time until a new measurement was done.
#' 
#' @param dataSemiLong Original data in semi-long format ("class data.frame").
#' @param timeColumn Character giving the column name of the observed times("logical vector"). It
#' is required that the observed times are discrete ("integer vector").
#' @param eventColumns Character vector giving the column names of the event
#' indicators (excluding censoring column)("character vector"). It is required that all events are
#' binary encoded. If the sum of all event indicators is zero, then this is
#' interpreted as a censored observation. Alternatively a column name of a
#' factor representing competing events can be given. In this case the argument
#' \emph{eventColumnsAsFactor} has to be set TRUE and the first level is assumed to
#' represent censoring.
#' @param eventColumnsAsFactor Should the argument eventColumns be intepreted
#' as column name of a factor variable ("logical vector")? Default is FALSE.
#' @param idColumn Name of column of identification number of persons as
#' character("character vector").
#' @param timeAsFactor Should the time intervals be coded as factor ("logical vector")? Default is
#' FALSE. In the default settings the discrete time intervals 
#' are treated as quantitative ("numeric vector").
#' @param responseAsFactor Should the response columns be given as factor ("logical vector")? 
#' Default is FALSE.
#' @return Original data set in long format with additional columns \itemize{
#' \item {obj} Gives identification number of objects (row index in short
#' format) (integer) \item {timeInt} Gives number of discrete time intervals
#' (factor) \item {responses} Columns with dimension count of events + 1
#' (censoring) \itemize{ \item {e0} No event (observation censored in specific
#' interval) \item {e1} Indicator of first event, 1 if event takes place and 0
#' otherwise \item ... ...  \item {ek} Indicator of last k-th event, 1 if event
#' takes place and 0 otherwise}
#' If argument responseAsFactor=TRUE, then responses will be coded as factor in one column.
#' }
#' @details In contrast to continuous survival (see e. g. \code{\link[survival]{Surv}}) 
#' the start and stop time notation is not used here. In discrete time survival analysis the only relevant
#' information is to use the stop time. Start time does not matter, because all discrete intervals need to be  
#' included in the long data set format to ensure consistent estimation. It is assumed that the supplied 
#' data set \emph{dataSemiLong} contains all repeated measurements of each cluster in semi-long format (e. g. persons). 
#' For further information see example \emph{Start-stop notation}.
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link{contToDisc}}, \code{\link{dataLong}},
#' \code{\link{dataLongCompRisks}}
#' @references
#' \insertRef{fahrmeirDiscSurv}{discSurv} \cr\cr
#' \insertRef{thompsonTreatment}{discSurv}
#' @keywords datagen
#' @examples
#' 
#' # Example Primary Biliary Cirrhosis data
#' library(survival)
#' pbcseq_example <- pbcseq
#' 
#' # Convert to months
#' pbcseq_example$day <- ceiling(pbcseq_example$day/30) + 1
#' names(pbcseq_example)[7] <- "month"
#' pbcseq_example$status <- factor(pbcseq_example$status)
#' 
#' # Convert to long format for time varying effects
#' pbcseq_exampleLong <- dataLongCompRisksTimeDep(dataSemiLong = pbcseq_example, timeColumn = "month", 
#' eventColumns = "status", eventColumnsAsFactor = TRUE, idColumn = "id", 
#' timeAsFactor = TRUE)
#' head(pbcseq_exampleLong)
#' 
#' #####################
#' # Start-stop notation
#' 
#' library(survival)
#' ?pbcseq
#' 
#' # Choose subset of patients
#' subsetID <- unique(pbcseq$id)[1:100]
#' pbcseq_mod <- pbcseq[pbcseq$id %in% subsetID, ]
#' 
#' # Convert to start stop notation
#' pbcseq_mod_split <- split(pbcseq_mod, pbcseq_mod$id)
#' pbcseq_mod_split <- lapply(1:length(pbcseq_mod_split), function(x) {
#' 
#'  cbind(pbcseq_mod_split[[x]], 
#'  start_time=c(0, pbcseq_mod_split[[x]][ - dim(pbcseq_mod_split[[x]])[1], "day"]),
#'  stop_time=pbcseq_mod_split[[x]][, "day"])
#'  
#' })
#' pbcseq_mod <- do.call(rbind, pbcseq_mod_split)
#' 
#' # Convert stop time to months
#' intervalDef <- c(quantile(pbcseq_mod$stop_time, probs = seq(0.1, 0.9, by=0.1)), Inf)
#' names(pbcseq_mod)
#' pbcseq_mod <- contToDisc(dataShort = pbcseq_mod, timeColumn = "stop_time", 
#'                          intervalLimits = intervalDef, equi = FALSE)
#' pbcseq_mod$status <- factor(pbcseq_mod$status)
#' 
#' # Conversion to data long format
#' pbcseq_mod_long <- dataLongCompRisksTimeDep(dataSemiLong = pbcseq_mod, timeColumn = "timeDisc", 
#'                                            eventColumns = "status",
#'                                           idColumn = "id", 
#'                                            eventColumnsAsFactor = TRUE, 
#'                                           responseAsFactor = TRUE,
#'                                           timeAsFactor = TRUE)
#' head(pbcseq_mod_long)
#' 
#' @export dataLongCompRisksTimeDep
dataLongCompRisksTimeDep <- function (dataSemiLong, timeColumn, eventColumns, 
                                      eventColumnsAsFactor = FALSE, idColumn, 
                                      timeAsFactor = FALSE, responseAsFactor = FALSE) {
  
  # Alternative input for competing risks coding 
  # as factor instead of multiple binary columns
  if(eventColumnsAsFactor) {
    respFact <- dataSemiLong[, eventColumns]
    respMat <- model.matrix(~., data.frame(respFact))[, -1]
    eventColumns <- levels(respFact)[-1]
    dimnames(respMat) [[2]] <- eventColumns
    dataSemiLong <- cbind(respMat, dataSemiLong)
  }
  
  # Construct censoring variable in short format
  indizeseventColumns <- sapply(1:length(eventColumns), function (x) which(eval(eventColumns [x]) == colnames(dataSemiLong)))
  NoeventColumns <- length(indizeseventColumns)
  # dataSet_eventColumns <- dataSemiLong [, indizeseventColumns]
  # eventColumnsorShort <- 1 - rowSums(dataSet_eventColumns)
  responsesShort <- cbind(1, matrix(0, nrow = nrow(dataSemiLong), ncol = NoeventColumns))
  namesEventColumnsNew <- paste("e", 0:NoeventColumns, sep="")
  dimnames(responsesShort) [[2]] <- namesEventColumnsNew
  dataSemiLong <- cbind(responsesShort, dataSemiLong)
  
  # Rearrange original data, that ID is in increasing order
  dataSemiLong <- dataSemiLong[order(dataSemiLong [, idColumn]), ]
  
  # Split data by persons
  splitDataSet <- split(x = dataSemiLong, f = dataSemiLong [, idColumn])
  lengthSplitDataSet <- 1:length(splitDataSet)
  
  # Get count of replicates in each split
  Counts <- lapply(splitDataSet, function (x) c(diff(x [, timeColumn]), 1))
  SumCounts <- as.numeric(sapply(Counts, function (x) sum(x)))
  
  # Replicate indices in each split
  Indizes <- lapply(lengthSplitDataSet, function (x) rep(x = 1:dim(splitDataSet [[x]])[1], times = Counts [[x]]))
  
  # Duplicate rows for each split and combine the results
  dataList <- lapply(lengthSplitDataSet, function (x) splitDataSet [[x]] [Indizes [[x]], ])
  dataList <- do.call(rbind, dataList)
  
  # Create ID variable in long format for each split and combine the results
  dataListID <- lapply(lengthSplitDataSet, function (x) rep(x, SumCounts [x]))
  dataListID <- do.call(c, dataListID)
  
  # Create time interval variable in long format for each split and combine results
  dataListTimeInt <- lapply(lengthSplitDataSet, function (x) 1:SumCounts [x])
  if(timeAsFactor) {
    dataListTimeInt <- factor(do.call(c, dataListTimeInt))
  }
  else{
    dataListTimeInt <- do.call(c, dataListTimeInt)
  }
  
  # Adjust responses
  cumSumSumCounts <- cumsum(SumCounts)
  rowSumsPre <- rowSums(dataList [, eventColumns])
  for(i in 1:length(cumSumSumCounts)) {
    if(rowSumsPre[ cumSumSumCounts[i] ] == 1) {
      # Replace censoring column with zero, if an event occurs
      dataList[cumSumSumCounts[i], namesEventColumnsNew[1] ] <- 0
      
      # Set respective event column to one
      relCol <- which(dataList [cumSumSumCounts[i], eventColumns] == 1) + 1
      dataList[cumSumSumCounts[i], namesEventColumnsNew[relCol] ] <- 1
    }
  }
  
  # Remove unnecessary columns
  if(eventColumnsAsFactor) {
    dataList <- dataList[, -which(names(dataList) %in% eventColumns)]
  }
  
  # Transform response columns to factor
  if(responseAsFactor) {
    responses <- dataList[, namesEventColumnsNew]
    colNameInd <- sapply(1:dim(responses)[1], 
                         function(x) which(responses[x, ] == 1) )
    responses <- factor(namesEventColumnsNew[colNameInd], levels = namesEventColumnsNew)
    dataList <- cbind(responses, dataList[, -c(1:length(namesEventColumnsNew))])
  }
  
  # Combine results and output
  dataComplete <- cbind(obj = dataListID, timeInt = dataListTimeInt, dataList)
  return(dataComplete)
}

######################################################################################
# Function to transform the univariate response in long format to Censoring variable
######################################################################################

# Description:
# Function for Transformation in censoring encoding in single event survival analysis
# NA values will be left out in glm! 

# Input
# dataSetLong: Original Data in long format
# respColumn: Variable name of discrete survival response as character
# timeColumn: Variable name of discrete time interval

# Output
# Original data.frame dataSetLong, but with added censoring process as first variable "yCens"

###############
# dataCensoring

#' Data Censoring Transformation for short formats
#' 
#' Function for transformation of discrete survival times in censoring
#' encoding. The original data is expanded to include the censoring process. 
#' Alternatively the long data format can also be augmented. With the new 
#' generated variable "yCens", the discrete censoring process can be analyzed 
#' instead of the discrete survival process. In discrete survival analysis this 
#' information is used to constructs weights for predictive evaluation measures. 
#' It is applicable in single event survival analysis.
#' 
#' @param dataShort Original data set in short format ("class data.frame").
#' @param eventColumns Name of event columns ("character vector"). The
#' event columns have to be in binary format. If the sum of all events equals
#' zero in a row, then this observation is interpreted as censored.
#' @param timeColumn Name of column with discrete time intervals ("character vector").
#' @param shortFormat Is the supplied data set \emph{dataShort} not preprocessed 
#' with function dataLong() ("logical vector")? Default is TRUE. If shortFormat=FALSE 
#' then it is assumed that the data set was augmented with function dataLong(). 
#' @return Original data set as argument \emph{dataShort}, but with added censoring
#' process as first variable in column "yCens".
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link{contToDisc}},
#' \code{\link{dataLong}}, \code{\link{dataLongTimeDep}},
#' \code{\link{dataLongCompRisks}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv} \cr\cr
#' \insertRef{fahrmeirDiscSurv}{discSurv} \cr\cr
#' \insertRef{thompsonTreatment}{discSurv}
#' @keywords datagen
#' @examples
#' 
#' library(pec)
#' data(cost)
#' head(cost)
#' IntBorders <- 1:ceiling(max(cost$time)/30)*30
#' subCost <- cost [1:100, ]
#' 
#' # Convert from days to months
#' CostMonths <- contToDisc(dataShort=subCost, timeColumn="time", intervalLimits=IntBorders)
#' head(CostMonths)
#' 
#' # Generate censoring process variable in short format
#' CostMonthsCensorShort <- dataCensoring (dataShort = CostMonths, 
#' eventColumns = "status", timeColumn = "time", shortFormat = TRUE)
#' head(CostMonthsCensorShort)
#' 
#' ################################
#' # Example with long data format
#' library(pec)
#' data(cost)
#' head(cost)
#' IntBorders <- 1:ceiling(max(cost$time)/30)*30
#' subCost <- cost [1:100, ]
#' 
#' # Convert from days to months
#' CostMonths <- contToDisc(dataShort = subCost, timeColumn = "time", intervalLimits = IntBorders)
#' head(CostMonths)
#' 
#' # Convert to long format based on months
#' CostMonthsLong <- dataLong(dataShort = CostMonths, timeColumn = "timeDisc", eventColumn = "status")
#' head(CostMonthsLong, 20)
#' 
#' # Generate censoring process variable
#' CostMonthsCensor <- dataCensoring (dataShort = CostMonthsLong, timeColumn = "timeInt", 
#' shortFormat = FALSE)
#' head(CostMonthsCensor)
#' tail(CostMonthsCensor [CostMonthsCensor$obj==1, ], 10)
#' tail(CostMonthsCensor [CostMonthsCensor$obj==3, ], 10)
#' 
#' @export dataCensoring
dataCensoring <- function(dataShort, eventColumns, timeColumn, shortFormat = TRUE){
  
  if(shortFormat){
    
    # Change relevant values of y to yCens
    ResultVec <- ifelse(rowSums(dataShort[, eventColumns, drop=FALSE]) == 0, 1, 0)
    
    # Exclude observations with time interval==1 and yCens==0
    deleteIndices <- ResultVec == 0 & dataShort[, timeColumn] == 1
    dataShort <- dataShort[!deleteIndices, ]
    ResultVec <- ResultVec[!deleteIndices]
    
    # Change time intervals for censoring process
    timeCens <- as.numeric(as.character(dataShort[, timeColumn]))
    timeCens[ResultVec == 0] <- timeCens[ResultVec == 0]-1
    
    # Append results to output
    NewDataSet <- cbind(yCens = ResultVec, timeCens = timeCens, dataShort)
    return(NewDataSet)
    
  } else{
    
    # Select last observed time
    lastTimeInd <- c((which(dataShort[, timeColumn] == 1)-1)[-1], dim(dataShort)[1])
    
    # Replicate original y
    ResultVec <- dataShort[, "y"]
    
    # Change relevant values of y to yCens
    ResultVec[lastTimeInd] <- ifelse(dataShort[lastTimeInd, "y"] == 0, 1, NA)
    
    # Append results to output
    NewDataSet <- cbind(yCens = ResultVec, dataShort)
    return(NewDataSet)
    
  }
}

#################
# dataLongSubdist

#' Data Matrix and Weights for Discrete Subdistribution Hazard Models
#' 
#' Generates the augmented data matrix and the weights required for discrete
#' subdistribution hazard modeling with right censoring.
#' 
#' This function sets up the augmented data matrix and the weights that are
#' needed for weighted maximum likelihood (ML) estimation of the discrete
#' subdistribution model proposed by Berger et al. (2018). The model is a
#' discrete-time extension of the original subdistribution model proposed by
#' Fine and Gray (1999).
#' 
#' @param dataShort Original data in short format ("class data.frame").
#' @param timeColumn Character specifying the column name of the observed event
#' times ("logical vector"). It is required that the observed times are discrete ("integer vector").
#' @param eventColumns Character vector specifying the column names of the
#' event indicators (excluding censoring events) ("logical vector"). It is required that a 0-1
#' coding is used for all events. The algorithm treats row sums of zero of all
#' event columns as censored.
#' @param eventColumnsAsFactor Should the argument \emph{eventColumns} be interpreted
#' as column name of a factor variable ("logical vector")? Default is FALSE.
#' @param eventFocus Column name of the event of interest (type 1 event) ("character vector").
#' @param timeAsFactor Logical indicating whether time should be coded as a
#' factor in the augmented data matrix("logical vector"). If FALSE, a numeric coding will be
#' used.
#' @param aggTimeFormat Instead of the usual long format, should every
#' observation have all time intervals? ("logical vector") Default is standard
#' long format. In the case of nonlinear risk score models, the time effect has
#' to be integrated out before these can be applied to the C-index.
#' @param lastTheoInt Gives the number of the last theoretic interval ("integer vector"). Only used, if aggTimeFormat==TRUE.
#' @return Data frame with additional column "subDistWeights". The latter
#' column contains the weights that are needed for fitting a weighted binary
#' regression model, as described in Berger et al. (2018). The weights are
#' calculated by a life table estimator for the censoring event.
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link{dataLong}}
#' @references 
#' \insertRef{bergerSubdist}{discSurv} \cr\cr
#' \insertRef{finePropHaz}{discSurv}
#' @keywords datagen
#' @examples
#' 
#' ################################
#' # Example with unemployment data
#' library(Ecdat)
#' data(UnempDur)
#' 
#' # Generate subsample, reduce number of intervals to k = 5
#' SubUnempDur <- UnempDur [1:500, ]
#' SubUnempDur$time <- as.numeric(cut(SubUnempDur$spell, c(0,4,8,16,28)))
#' 
#' # Convert competing risks data to long format
#' # The event of interest is re-employment at full job
#' SubUnempDurLong <- dataLongSubDist (dataShort=SubUnempDur, timeColumn = "time", 
#' eventColumns=c("censor1", "censor2", "censor3"), eventFocus="censor1")
#' head(SubUnempDurLong)
#' 
#' # Fit discrete subdistribution hazard model with logistic link function
#' logisticSubDistr <- glm(y ~ timeInt + ui + age + logwage,
#'                     family=binomial(), data = SubUnempDurLong, 
#'                     weights = SubUnempDurLong$subDistWeights)
#' summary(logisticSubDistr)
#' 
#' ########################################
#' # Simulation 
#' # Discrete subdistribution hazards model
#' 
#' # Simulate covariates as multivariate normal distribution
#' library(mvnfast)
#' set.seed(1980)
#' X <- mvnfast::rmvn(n = 1000, mu = rep(0, 4), sigma = diag(4))
#'
#' # Specification of two discrete cause specific hazards with four intervals
#' # Event 1
#' theoInterval <- 4
#' betaCoef_event1 <- seq(-1, 1, length.out = 5)[-3]
#' timeInt_event1 <- seq(0.1, -0.1, length.out = theoInterval-1)
#' linPred_event1 <- c(X %*% betaCoef_event1)
#' # Event 2
#' betaCoef_event2 <- seq(-0.5, 0.5, length.out = 5)[-3]
#' timeInt_event2 <- seq(-0.1, 0.1, length.out = theoInterval-1)
#' linPred_event2 <- c(X %*% betaCoef_event2)
#' # Discrete cause specific hazards in last theoretical interval
#' theoHaz_event1 <- 0.5
#' theoHaz_event2 <- 0.5
#'
#' # Derive discrete all cause hazard
#' haz_event1_X <- cbind(sapply(1:length(timeInt_event1), 
#'                             function(x) exp(linPred_event1 + timeInt_event1[x]) / 
#'                               (1 + exp(linPred_event1 + timeInt_event1[x]) + 
#'                               exp(linPred_event2 + timeInt_event2[x])) ), 
#'                      theoHaz_event1)
#'
#' haz_event2_X <- cbind(sapply(1:length(timeInt_event2), 
#'                             function(x) exp(linPred_event2 + timeInt_event2[x]) / 
#'                               (1 + exp(linPred_event1 + timeInt_event1[x]) + 
#'                               exp(linPred_event2 + timeInt_event2[x]) ) ),
#'                      theoHaz_event2)
#' allCauseHaz_X <- haz_event1_X + haz_event2_X
#'
#' # Derive discrete cumulative incidence function of event 1 given covariates
#' p_T_event1_X <- haz_event1_X * cbind(1, (1-allCauseHaz_X)[, -dim(allCauseHaz_X)[2]])
#' cumInc_event1_X <-  t(sapply(1:dim(p_T_event1_X)[1], function(x) cumsum(p_T_event1_X[x, ])))
#'
#' # Calculate all cause probability P(T=t | X)
#' pT_X <- t(sapply(1:dim(allCauseHaz_X)[1], function(i) estMargProb(allCauseHaz_X[i, ]) ))
#'
#' # Calculate event probability given time interval P(R=r | T=t, X)
#' pR_T_X_event1 <- haz_event1_X / (haz_event1_X + haz_event2_X)
#'
#' # Simulate discrete survival times
#' survT <- sapply(1:dim(pT_X)[1], function(i) sample(x = 1:(length(timeInt_event1)+1), 
#'                                                    size = 1, prob = pT_X[i, ]) )
#' censT <- sample(x = 1:(length(timeInt_event1)+1), size = dim(pT_X)[1], 
#'                prob = rep(1/(length(timeInt_event1) + 1), (length(timeInt_event1) + 1)), 
#'                replace = TRUE)
#'
#' # Calculate observed times
#' obsT <- ifelse(survT <= censT, survT, censT)
#' obsEvent <- rep(0, length(obsT))
#' obsEvent <- sapply(1:length(obsT), 
#'                   function(i) if(survT[i] <= censT[i]){
#'                     return(sample(x = c(1, 2), size = 1, 
#'                     prob = c(pR_T_X_event1[i, obsT[i]  ], 
#'                     1 - pR_T_X_event1[i, obsT[i]  ]) ))
#'                   } else{
#'                     
#'                     return(0)
#'                   }
#')
#'
#' # Recode last interval to censored
#' lastInterval <- obsT == theoInterval
#' obsT[lastInterval] <- theoInterval-1
#' obsEvent[lastInterval] <- 0
#' obsT <- factor(obsT)
#' obsEvent <- factor(obsEvent)
#'
#' # Data preparation
#' datShort <- data.frame(event = factor(obsEvent), time=obsT, X)
#'
#' # Conversion to long data format
#' datLongSub <- dataLongSubDist(dataShort = datShort, timeColumn = "time",
#'                              eventColumns = "event", eventFocus = 1, eventColumnsAsFactor = TRUE)
#'
#' # Estimate discrete subdistribution hazard model
#' estSubModel <- glm(formula = y ~ timeInt + X1 + X2 + X3 + X4, data = datLongSub,
#'                   family = binomial(link = "logit"), weights = datLongSub$subDistWeights)
#'
#' # Predict cumulative incidence function of first event
#' predSubHaz1 <- predict(estSubModel, newdata = datLongSub[datLongSub$obj == 2, ], type = "response")
#' mean(((1 - estSurv(predSubHaz1)) - cumInc_event1_X[2, 1:3])^2)
#' 
#' @export dataLongSubDist
dataLongSubDist <- function(dataShort, timeColumn, eventColumns, 
                             eventColumnsAsFactor=FALSE, eventFocus, timeAsFactor=FALSE,
                            aggTimeFormat=FALSE, lastTheoInt=NULL) {
  #######################
  # Construct long format
  
  # Alternative input for competing risks coding 
  # as factor instead of multiple binary columns
  if(eventColumnsAsFactor) {
    respFact <- dataShort[, eventColumns]
    respMat <- model.matrix(~., data.frame(respFact))[, -1]
    eventColumns <- levels(respFact)[-1]
    dimnames(respMat) [[2]] <- eventColumns
    dataShort <- cbind(respMat, dataShort)
  }
  
  # Reorder event columns
  # First column correspond to the event of focus
  eventColumns <- c(eventColumns[eventColumns == eventFocus], 
                    eventColumns[eventColumns != eventFocus])

  # dataLongCompRisks
  # # Construct object counter, covariates in long format and discrete time intervals
  # if(aggTimeFormat){
  #   
  #   obj <- rep(1:nrow(dataShort), each = lastTheoInt)
  #   dataSetLong <- dataShort[obj,]
  #   timeInt <- rep(1:lastTheoInt, nrow(dataShort))
  #   if(timeAsFactor){timeInt <- factor(timeInt)}
  #   
  # }
  
  # Construct object counter
  dataSet_timeColumn_org <- as.numeric(as.character(dataShort[, timeColumn]))
  
  if(aggTimeFormat){
    tmax <- lastTheoInt
  } else{
    tmax <- max(dataSet_timeColumn_org)
  }

  delta_i <- rowSums(dataShort[, eventColumns]) == 1
  epsilon_i <- sapply(1:dim(dataShort)[1], 
                      function(x) ifelse(delta_i[x], 
                                         which(dataShort[x, eventColumns] == 1), 0))
  dataSet_timeColumn <- ifelse(delta_i & 
                                 dataShort[, eventFocus] == 0, tmax, 
                               dataSet_timeColumn_org)
  obj <- rep(1:nrow(dataShort), each = tmax)
  
  # Construct time intervals
  timeInt <- c(sapply(1:length(dataSet_timeColumn), 
                      function(k) 1:tmax))
  
  # dataLongCompRisks
  # # Construct responses and censoring variable
  # if(aggTimeFormat){
  #   
  #   NoeventColumns <- length(indizeseventColumns)
  #   responses <- matrix(0, nrow = 0, ncol = NoeventColumns + 1)
  #   for(i in 1:dim(dataShort) [1]) {
  #     row_rep <- as.numeric(as.character(dataShort [i, indextimeColumn])) - 1
  #     mat_temp <- matrix(rep(c(1, rep(0, NoeventColumns)), row_rep), nrow = row_rep, ncol = NoeventColumns + 1, byrow = TRUE)
  #     mat_temp <- rbind(mat_temp, responsesShort [i, ])
  #     mat_temp_after <- matrix(rep(c(1, rep(0, NoeventColumns)), lastTheoInt - row_rep - 1), 
  #                              nrow = lastTheoInt - row_rep - 1, ncol = NoeventColumns + 1, byrow = TRUE)
  #     responses <- rbind(responses, mat_temp, mat_temp_after)
  #   }
  #   dimnames(responses) [[2]] <- paste("e", 0:NoeventColumns, sep = "")
  #   
  # }
  
  # Construct response for event of interest
  eventFocus <- dataShort[, eventFocus]
  y <- c(sapply(1:length(dataSet_timeColumn), 
                function(k) c(rep(0, dataSet_timeColumn[k]-1), 
                              eventFocus[k], 
                              rep(0, tmax-dataSet_timeColumn[k]) ) ))
  
  # Estimation of weights
  estG <- estSurvCens(dataShort = dataShort, timeColumn = timeColumn, 
                      eventColumns = eventColumns)
  weights <- estG[timeInt] / estG[pmin(dataSet_timeColumn_org[obj], timeInt)] * 
    (ifelse(timeInt <= dataSet_timeColumn_org[obj], 1, 0) + 
       ifelse(dataSet_timeColumn_org[obj] <= (timeInt - 1) & 
                delta_i[obj] * epsilon_i[obj] > 1, 1, 0))
  
  # Combine results
  if(timeAsFactor){timeInt <- factor(timeInt)}
  dataSetLong <- cbind(obj = obj, timeInt = timeInt, y = y, dataShort[obj, ])
  Output <- cbind(dataSetLong, subDistWeights = weights)
  return(Output)
} 

#################################
# Multi Spell data transformation

#' Data long transformation for multi spell analysis
#' 
#' Transform data from short format into long format for discrete multi spell
#' survival analysis and right censoring.
#' 
#' If the data has continuous survival times, the response may be transformed
#' to discrete intervals using function \code{\link{contToDisc}}. The discrete
#' time variable needs to be strictly increasing for each person, because
#' otherwise the order of the events is not distinguishable. Here is an example
#' data structure in short format prior augmentation with three possible
#' states: \ idColumn=1, 1, ... , 1, 2, 2, ... , n \ timeColumn= t_ID1_1 <
#' t_ID1_1 < ... < t_ID1_k, t_ID2_1 < t_ID2_2 < ... < t_ID2_k, ... \
#' eventColumn = 0, 1, ... , 2, 1, 0, ... , 0
#' 
#' @param dataSemiLong Original data in semi-long format ("class data.frame").
#' @param timeColumn Character giving the column name of the observed times. It
#' is required that the observed times are discrete ("character vector").
#' @param eventColumn Column name of the event status ("character vector"). The
#' events can take multiple values on a discrete scale (0, 1, 2, ...) and
#' repetition of events is allowed (integer vector or class factor). 
#' It is assumed that the number zero corresponds to censoring and all number > 0 
#' represent the observed states between transitions.
#' @param idColumn Name of column of identification number of persons as
#' character("character vector").
#' @param timeAsFactor Should the time intervals be coded as factor ("logical vector")? Default is FALSE. 
#' In the default settings the discrete time intervals are treated as quantitative ("numeric vector").
#' @param spellAsFactor Should the spells be coded as factor ("logical vector")? Default is
#' not to use factor. If the argument is false, the column is coded as numeric.
#' @return Original data.frame with three additional columns: \itemize{ \item
#' {obj} Index of persons as integer vector \item {timeInt} Index of time
#' intervals (factor or integer vector) \item {spell} The spell gives the actual 
#' state of each individual within a given discrete interval.
#' \item {e0} Response transition in long format as binary vector. Column \emph{e0} represents censoring. 
#' If \emph{e0} is coded one in the in the last observed time interval \emph{timeInt} of a person, 
#' then this observation was censored. \item {e1} Response in long format as binary vector. 
#' The column \emph{e1} represents the transition to the first event state. 
#' \item {eX} Response in long format as binary vector. 
#' The column \emph{eX} represents the transition to the last event state out of the set of possible states 
#' "1, 2, 3, ..., X". \item ... Expanded columns of original data set.}
#' @details The starting state of each individual is assumed to given with time interval
#' equals zero. For example in an illness-death model with three states ("healthy", "illness", "death")
#' if an individual was healthy at the beginning of the study this has to be encoded with
#' discrete time interval set to zero and event state "healthy".
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' @seealso \code{\link{contToDisc}}, \code{\link{dataLongTimeDep}},
#' \code{\link{dataLongCompRisks}}, \code{\link{dataLongCompRisks}}
#' @references 
#' \insertRef{tutzModelDisc}{discSurv} \cr\cr
#' \insertRef{fahrmeirDiscSurv}{discSurv} \cr\cr
#' \insertRef{thompsonTreatment}{discSurv}
#' @keywords datagen
#' @examples
#' 
#' ################################
#' # Example with unemployment data
#' data(unempMultiSpell)
#' 
#' # Select subsample of first 500 persons
#' unempSub <- unempMultiSpell[unempMultiSpell$id %in% 1:250,]
#' 
#' # Expansion from semi-long to long format
#' unempLong <- dataLongMultiSpell(dataSemiLong=unempSub, timeColumn = "year",
#'                                 eventColumn="spell", idColumn="id", 
#'                                 spellAsFactor=TRUE, timeAsFactor=FALSE)
#'
#' head(unempLong, 25)
#'
#' # Fit discrete multi-state model regression model
#' library(VGAM)
#'
#' model <- vgam(cbind(e0, e1, e2, e3, e4) ~ 0 + s(timeInt) + age:spell, 
#' data = unempLong, family = multinomial(refLevel="e0"))
#'              
#' ############################
#' # Example with artificial data
#' 
#' # Seed specification
#' set.seed(-2578)
#' 
#' # Construction of data set
#' # Censoring and three possible states (0, 1, 2, 3)
#' # Discrete time intervals (1, 2, ... , 10)
#' # Noninfluential variable x ~ N(0, 1)
#' datFrame <- data.frame(
#'  ID = c(rep(1, 6), rep(2, 4), rep(3, 3), rep(4, 2), rep(5, 4), 
#'       rep(6, 5), rep(7, 7), rep(8, 8)),
#'  time = c(c(0, 2, 5, 6, 8, 10), c(0, 1, 6, 7), c(0, 9, 10), c(0, 6), c(0, 2, 3, 4), 
#'         c(0, 3, 4, 7, 9), c(0, 2, 3, 5, 7, 8, 10), c(0, 1, 3, 4, 6, 7, 8, 9) ),
#'  state = c(c(2, 1, 3, 2, 1, 0), c(3, 1, 2, 2), c(2, 2, 1), c(1, 2), c(3, 2, 2, 0), 
#'          c(1, 3, 2, 1, 3), c(1, 1, 2, 3, 2, 1, 3), c(3, 2, 3, 2, 1, 1, 2, 3) ),
#'  x = rnorm(n=6+4+3+2+4+5+7+8) )
#'
#' # Transformation to long format
#' datFrameLong <- dataLongMultiSpell(dataSemiLong=datFrame, timeColumn="time",
#'                                    eventColumn="state", idColumn="ID", 
#'                                    spellAsFactor=TRUE)
#' head(datFrameLong, 25)

# Fit multi state model without autoregressive terms
#' library(VGAM)
#' cRm <- vglm(cbind(e0, e1, e2, e3) ~ 0 + timeInt + x:spell, 
#' data = datFrameLong, family = "multinomial")
#' summary(cRm)
# -> As expected there is no significant effect of x across all different spells.
#' 
#' @export dataLongMultiSpell
dataLongMultiSpell <- function(dataSemiLong, timeColumn, eventColumn, idColumn,
                               timeAsFactor=FALSE, spellAsFactor=FALSE){
  
  # Add column starting spell (constant on person-level)
  dataSemiLong$startSpell <- rep(as.numeric(as.character(
    dataSemiLong[dataSemiLong[, timeColumn] == 0, eventColumn])), 
    times=table(dataSemiLong[, idColumn]))
  dataSemiLong <- dataSemiLong[dataSemiLong[,timeColumn] != 0,]
  
  # Expand responses
  # Assumption: Reference level corresponds to censoring
  if(!is.factor(dataSemiLong[, eventColumn])){
    dataSemiLong[, eventColumn] <- factor(dataSemiLong[, eventColumn])
  }
  eventExpand <- model.matrix(~.-1, dataSemiLong[, eventColumn, drop = FALSE])
  dimnames(eventExpand)[[2]] <- paste("e", 0:(dim(eventExpand)[2]-1), sep = "")
  dataSemiLong <- cbind(obj = NA, timeInt = NA, spell = dataSemiLong[, "startSpell"], 
                          eventExpand, dataSemiLong)
  
  # Extract indices for IDs from time
  splitDat <- split(dataSemiLong, dataSemiLong[, idColumn])
  lenSplitDat <- length(splitDat)
  
  # Find discrete time intervals with time gaps between spells
  for( j in 1:lenSplitDat ){
    
    # Calculate time differences
    actualTime <- splitDat[[j]][, timeColumn]
    diffTime <- diff( actualTime )
    # Define standard indices
    defIndices <- 1:length(actualTime)
    if(actualTime[1] > 1){
      defIndices <- c( defIndices, rep(1, actualTime[1] - 1) )
    }
    
    for( i in 1:length(diffTime) ){
        if(length(diffTime) > 0){
          if(diffTime[i] > 1){
            defIndices <- c( defIndices, rep(i+1, diffTime[i] - 1) )
          }
        }
    }
    # Enlarge data set
    defIndices <- sort(defIndices)
    splitDat[[j]] <- splitDat[[j]][defIndices, ]
    splitDat[[j]][, "obj"] <- j
    splitDat[[j]][, "timeInt"] <- 1:dim(splitDat[[j]])[1]

    # Adapt responses
    for( k in 1:length(actualTime) ){
      
      if(length(splitDat[[j]][defIndices == k, "e0"]) > 1){
        relInd <- which(defIndices == k)
        relInd <- relInd[-length(relInd)]
        splitDat[[j]][relInd, "e0"] <- 1
        splitDat[[j]][relInd, paste("e", 
                                    1:(dim(eventExpand)[2]-1), sep = "")] <- 0
        
        
      }
      
    }
    
    # Add current spell variable
    rowSumsEvent <- rowSums(splitDat[[j]][, paste("e", 1:(dim(eventExpand)[2]-1), sep = "")])
    rowSumsIndices <- which(rowSumsEvent == 1)
    maxObsTime <- max(dataSemiLong[, timeColumn])
    if( tail(splitDat[[j]][, "e0"], 1) == 1 ){
      
      rowSumsIndices <- c(rowSumsIndices, maxObsTime)
      
    }
    indicesSpell <- lapply(0:(length(rowSumsIndices)-1), 
                           function(x) {
                             
                             if(x > 0){
                               
                               return((rowSumsIndices[x] + 1):rowSumsIndices[x + 1])
                               
                             } else{
                               
                               return(1:rowSumsIndices[1])
                               
                             }})
    
    if( length(indicesSpell) > 1 ){
      
      for( l in 2:length(indicesSpell) ){
        
        splitDat[[j]][ indicesSpell[[l]], "spell" ] <- as.numeric(as.character(
          splitDat[[j]][ indicesSpell[[l - 1]], eventColumn][1]))
        
      }
      
    }

  }
  
  # Return full data set
  allDat <- as.data.frame(rbindlist(splitDat))
  
  # Should time be formatted as factor?
  if(timeAsFactor){
    allDat$timeInt <- factor(allDat$timeInt)
  }
  if(spellAsFactor){
    allDat$spell <- factor(allDat$spell)
  }
  
  # Data clearing
  allDat <- allDat[!rowSums(is.na(allDat[, 1:2]))>0,]
  allDat <- allDat[, -which(names(allDat)=="startSpell")]
  
  return(allDat)
}
