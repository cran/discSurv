#########################################################
# Transform continuous time variables to a discrete grid
#########################################################

# Description
# Discretizes continuous time variable into a specified grid of correct censored data

# Input
# dataSet: Data frame with observed times, event indicator and covariates
# timeColumn: Character giving the column name of the observed times
# intervalLimits: Numeric vector of the right interval borders, e. g. if the intervals are
  # [0, a_1), [a_1, a_2), [a_2, a_{\max}), then intervalLimits = c(a_1, a_2, a_{\max})

# Output
# Original data frame with discretized response in the first column

contToDisc <- function(dataSet, timeColumn, intervalLimits, equi=FALSE) {
  
  # Input checks
  # Correct formats
  if(!is.data.frame(dataSet)) {stop("Argument *dataSet* is not in the correct format! Please specify as data.frame object.")}
  if(!(is.character(timeColumn))) {stop("Argument *timeColumn* is not in the correct format! Please specify as character.")}
  if(length(timeColumn)!=1) {stop("Argument *timeColumn* is not in the correct format! Please specify as scalar with length one.")}
  if(equi==FALSE) {
    if(!(is.vector(intervalLimits))) {stop("Argument *intervalLimits* is not in the correct format! Please specify as numeric vector.")}
    if(!(is.numeric(intervalLimits))) {stop("Argument *intervalLimits* is not in the correct format! Please specify as numeric vector.")}
  }
  if(equi==TRUE & length(intervalLimits) != 1) {stop("Argument *intervalLimits* is not in the correct format! Please specify an scalar integer value.")}
  if(equi==TRUE & any(floor(intervalLimits) != intervalLimits)) {stop("Argument *intervalLimits* is not in the correct format! Please specify an scalar integer value.")}
  
  # Can *timeColumn* be accessed in *dataSet*?
  if(any(class(tryCatch(dataSet [, timeColumn], error=function (e) e)) == "error")) 
  {stop("*timeColumn* is not available in *dataSet*! Please specify the correct column of observed times.")}
  
  #####################
  # Main Code
  
  # Extract observed time
  obsTime <- dataSet [, timeColumn]
  
  if(equi==FALSE) {
    # Include zero in interval limits
    intervalLimits <- c(0, intervalLimits)
    
    # Check if breaks are not unique and add little bit of jitter in that case
    # Then the breaks are unique!
    Check <- tryCatch(cut(obsTime, intervalLimits, right = FALSE), error=function (e) "Bug!")
    if(identical(Check, "Bug!")) {
      intervalLimits [-1] <- jitter(intervalLimits [-1], factor=0.01)
    }
  }
  
  # Divide the time to discrete intervals
  timeDisc <- as.ordered(as.numeric(cut(obsTime, intervalLimits, right = FALSE)))
  
  # Include discretized time in in original data set
  dataSet <- cbind(timeDisc=timeDisc, dataSet)
  return(dataSet)
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
# dataSet: Data frame with observed times, event indicator and covariates
# timeColumn: Character giving the column name of the observed times. 
  # It is required that the observed times are discrete (integer)
# censColumn: Character giving the column name of the event indicator
 # It is required that this is a binary variable with 1=="event" and 0=="censored"

# Output: Original data.frame with three additional columns
# obj: Index of persons as integer vector
# timeInt: Index of time intervals as integer vector
# y: Response in long format as binary vector. 1=="event happens in period timeInt" and 0 otherwise

dataLong <- function (dataSet, timeColumn, censColumn, timeAsFactor=TRUE,
                      remLastInt=FALSE, aggTimeFormat=FALSE, lastTheoInt=NULL) {

  # Input checks
  if(!is.data.frame(dataSet)) {stop("Argument *dataSet* is not in the correct format! Please specify as data.frame object.")}
  if(!(is.character(timeColumn))) {stop("Argument *timeColumn* is not in the correct format! Please specify as character.")}
  if(length(timeColumn)!=1) {stop("Argument *timeColumn* is not in the correct format! Please specify as scalar with length one.")}
  if(!(is.character(censColumn))) {stop("Argument *censColumn* is not in the correct format! Please specify as character.")}
  if(length(censColumn)!=1) {stop("Argument *censColumn* is not in the correct format! Please specify as scalar with length one.")}

  # Can *timeColumn* be accessed in *dataSet*?
  if(any(class(tryCatch(dataSet [, timeColumn], error=function (e) e)) == "error")) {
    stop("*timeColumn* is not available in *dataSet*! Please specify the correct column of observed times.")
  }
  
  # Can *censColumn* be accessed in *dataSet*?
  if(any(class(tryCatch(dataSet [, censColumn], error=function (e) e)) == "error")) {
    stop("*censColumn* is not available in *dataSet*! Please specify the correct column of the event indicator.")
  }
  
  # Check if observed times are only integer values
  if(!all(dataSet [, timeColumn]==floor(as.numeric(as.character(dataSet [, timeColumn]))))) {
    stop("*timeColumn* has not only integer values! Please convert the observed time in discrete intervals.")
  }
  
  # Check if event indicator has only zero or one values
  if(!all(as.numeric(as.character(dataSet [, censColumn]))==0 | as.numeric(as.character(dataSet [, censColumn]))==1)) {
    stop("*censColumn* is not a binary vector! Please check, that events equals 1 and 0 otherwise.")
  }
  
  ##################
  # Main Code
  
  # Extract column index of survival and censoring times
  c1 <- which(eval(timeColumn)==colnames(dataSet))
  c2 <- which(eval(censColumn)==colnames(dataSet))
  
  # Construct indices of persons
  if(aggTimeFormat){
    obj <- rep(1:nrow(dataSet), each=lastTheoInt)
  } else{
    obj <- rep(1:nrow(dataSet), as.vector(dataSet[, c1]))
  }

  # Long format of covariates
  dataSetLong <- dataSet[obj, ]
  
  # Calculate discrete time interval
  if(aggTimeFormat){
    if(timeAsFactor) {
      timeInt <- factor( rep(1:lastTheoInt, nrow(dataSet)) )
    }
    else{
      timeInt <- rep(1:lastTheoInt, nrow(dataSet))
    }
    # Calculate response
    y <- c(unlist(apply(dataSet, 1, FUN = 
        function(k) c(rep(0, as.numeric(k[c1])-1), as.numeric(k[c2]), 
                      rep(0, lastTheoInt - as.numeric(k[c1]) ) ) )))
    
  } else{
    if(timeAsFactor) {
      timeInt <- factor(unlist(apply(dataSet, 1, FUN = function(k) 1:k[c1])))
    }
    else{
      timeInt <- unlist(apply(dataSet, 1, FUN = function(k) 1:k[c1]))
    }
    # Calculate response
    y <- c(unlist(apply(dataSet, 1, FUN = 
         function(k) c(rep(0, as.numeric(k[c1])-1), as.numeric(k[c2])) )))
  }
  
  # Aggregate results in one data.frame
  dataSetLong <- cbind(obj, timeInt, y, dataSetLong)
  
  if(remLastInt) {
    # Remove cases with observed values in last interval, 
    # because the hazard is always 1
    # t==t_max & event==1 -> t==t_max-1 & event==0
    remInd <- which(dataSetLong$y==1 & 
                      dataSetLong$timeInt==max(as.numeric(as.character(
                        dataSetLong$timeInt))) )
    if(length(remInd)!=0) {
      dataSetLong <- dataSetLong[-remInd, ]
    }
  }
  
  return(dataSetLong)
} 

#########################################################
# Datalong for time dependent covariates
#########################################################

# Description:
# Transform data in short format into long format for discrete survival analysis
# Data is assumed to include time varying covariates, e. g. follow up visits are allowed
# It is assumed that the covariates stay constant between measurements, where no information is available

# Input
# dataSet: Data frame with observed times, event indicator and covariates
# timeColumn: Character giving the column name of the observed times. 
# It is required that the observed times are discrete.
# censColumn: Character giving the column name of the event indicator
# It is required that this is a binary variable with 1=="event" and 0=="censored"
# idColumn: Gives the identification number of the persons

# Output: Original data.frame with three additional columns
# obj: Index of persons as integer vector
# timeInt: Index of time intervals as integer vector
# y: Response in long format as binary vector. 1=="event happens in period timeInt" and 0 otherwise

dataLongTimeDep <- function (dataSet, timeColumn, censColumn, idColumn, timeAsFactor=TRUE) {

  # Input checks
  if(!is.data.frame(dataSet)) {stop("Argument *dataSet* is not in the correct format! Please specify as data.frame object.")}
  if(!(is.character(timeColumn))) {stop("Argument *timeColumn* is not in the correct format! Please specify as character.")}
  if(length(timeColumn)!=1) {stop("Argument *timeColumn* is not in the correct format! Please specify as scalar with length one.")}
  if(!(is.character(censColumn))) {stop("Argument *censColumn* is not in the correct format! Please specify as character.")}
  if(length(censColumn)!=1) {stop("Argument *censColumn* is not in the correct format! Please specify as scalar with length one.")}
  
  # Can *timeColumn* be accessed in *dataSet*?
  if(any(class(tryCatch(dataSet [, timeColumn], error=function (e) e)) == "error")) {
    stop("*timeColumn* is not available in *dataSet*! Please specify the correct column of observed times.")
  }
  # Can *censColumn* be accessed in *dataSet*?
  if(any(class(tryCatch(dataSet [, censColumn], error=function (e) e)) == "error")) {
    stop("*censColumn* is not available in *dataSet*! Please specify the correct column of the event indicator.")
  }
  # Can *idColumn* be accessed in *dataSet*?
  if(any(class(tryCatch(dataSet [, idColumn], error=function (e) e)) == "error")) {
    stop("*idColumn* is not available in *dataSet*! Please specify the correct column of the identification number.")
  }
  
  # Check if observed times are only integer values
  if(!all(dataSet [, timeColumn]==floor(as.numeric(as.character(dataSet [, timeColumn]))))) {
    stop("*timeColumn* has not only integer values! Please convert the observed time in discrete intervals.")
  }
  # Check if event indicator has only zero or one values
  if(!all(as.numeric(as.character(dataSet [, censColumn]))==0 | as.numeric(as.character(dataSet [, censColumn]))==1)) {
    stop("*censColumn* is not a binary vector! Please check, that events equals 1 and 0 otherwise.")
  }
  
  ###############
  # Main code

  # Rearrange original data, that ID is in increasing order
  dataSet <- dataSet[order(dataSet [, idColumn]), ]
  
  # Split data by persons
  splitDataSet <- split(x=dataSet, f=dataSet [, idColumn])
  lengthSplitDataSet <- 1:length(splitDataSet)
  
  # Get count of replicates in each split
  Counts <- lapply(splitDataSet, function (x) c(diff(x [, timeColumn]), 1))
  SumCounts <- as.numeric(sapply(Counts, function (x) sum(x)))
  
  # Replicate indices in each split
  Indizes <- lapply(lengthSplitDataSet, function (x) rep(x=1:dim(splitDataSet [[x]])[1], times=Counts [[x]]))
  
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
  dataListResponse <- lapply(lengthSplitDataSet, function (x) c(rep(0, SumCounts [x]-1), as.numeric(as.character(tail(splitDataSet [[x]] [, censColumn], 1)))))
  dataListResponse <- do.call(c, dataListResponse)
  
  # Combine overall results and output
  dataComplete <- cbind(obj=dataListID, timeInt=dataListTimeInt, y=dataListResponse, dataList)
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
# dataSet: Complete data set in short format (data.frame)
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
# dataSet: Original data with replicated rows

dataLongCompRisks <- function (dataSet, timeColumn, eventColumns, 
                               eventColumnsAsFactor=FALSE, timeAsFactor=TRUE) {

  # Input checks
  if(!is.data.frame(dataSet)) {stop("Argument *dataSet* is not in the correct format! Please specify as data.frame object.")}
  if(!(is.character(timeColumn))) {stop("Argument *timeColumn* is not in the correct format! Please specify as character.")}
  if(length(timeColumn)!=1) {stop("Argument *timeColumn* is not in the correct format! Please specify as scalar with length one.")}
  if(!all(is.character(eventColumns))) {stop("Argument *eventColumns* is not in the correct format! Please specify as character.")}
  
  # Can *timeColumn* be accessed in *dataSet*?
  if(any(class(tryCatch(dataSet [, timeColumn], error=function (e) e)) == "error")) {
    stop("*timeColumn* is not available in *dataSet*! Please specify the correct column of observed times.")
  }
  # Can *eventColumns* be accessed in *dataSet*?
  if(any(class(tryCatch(dataSet [, eventColumns], error=function (e) e)) == "error")) {
    stop("*eventColumns* is not available in *dataSet*! Please specify the correct column of the event indicator.")
  }
  
  # Check if observed times are only integer values
  if(!all(dataSet [, timeColumn]==floor(as.numeric(as.character(dataSet [, timeColumn]))))) {
    stop("*timeColumn* has not only integer values! Please convert the observed time in discrete intervals.")
  }

  #############
  # Main Code
  
  # Alternative input for competing risks coding 
  # as factor instead of multiple binary columns
  if(eventColumnsAsFactor) {
    respFact <- dataSet[, eventColumns]
    respMat <- model.matrix(~., data.frame(respFact))[, -1]
    eventColumns <- levels(respFact)[-1]
    dimnames(respMat) [[2]] <- eventColumns
    dataSet <- cbind(respMat, dataSet)
  }
  
  # Check if event indicators have only zero or one values
  checkVec <- vector("logical", length(eventColumns))
  for(i in 1:length(eventColumns)) {
    checkVec [i] <- all(as.numeric(as.character(dataSet [, eventColumns [i] ]))==0 | as.numeric(as.character(dataSet [, eventColumns [i] ]))==1)
  }
  if(!all(checkVec)) {
    stop("*eventColumns* is not a binary vector! Please check, that events equals 1 and 0 otherwise.")
  }
  
  # Index of columns of timeColumn and event indicators
  indextimeColumn <- which(eval(timeColumn)==colnames(dataSet))
  indizeseventColumns <- sapply(1:length(eventColumns), function (x) which(eval(eventColumns [x])==colnames(dataSet)))
  dataSet_timeColumn <- as.numeric(as.character(dataSet[, indextimeColumn]))
  
  # Construct object counter, covariates in long format and discrete time intervals
  obj <- rep(1:nrow(dataSet), dataSet_timeColumn)
  dataSetLong <- dataSet[obj,]
  timeInt <- unlist(apply(dataSet, 1, FUN = function(k) 1:k[indextimeColumn]))
  if(timeAsFactor){timeInt <- factor(timeInt)}
  
  # Construct censoring variable in short format
  dataSet_eventColumns <- dataSet [, indizeseventColumns]
  eventColumnsorShort <- 1 - rowSums(dataSet_eventColumns)
  responsesShort <- as.matrix(cbind(eventColumnsorShort, dataSet_eventColumns))
  dimnames(responsesShort) [[2]] <- 1:length(dimnames(responsesShort) [[2]])
  
  # Construct responses and censoring variable
  NoeventColumns <- length(indizeseventColumns)
  responses <- matrix(0, nrow=0, ncol=NoeventColumns + 1)
  for(i in 1:dim(dataSet) [1]) {
    row_rep <- as.numeric(as.character(dataSet [i, indextimeColumn])) - 1
    mat_temp <- matrix(rep(c(1, rep(0, NoeventColumns)), row_rep), nrow=row_rep, ncol=NoeventColumns + 1, byrow=TRUE)
    mat_temp <- rbind(mat_temp, responsesShort [i, ])
    responses <- rbind(responses, mat_temp)
  }
  dimnames(responses) [[2]] <- paste("e", 0:NoeventColumns, sep="")
  
  # Remove unnecessary columns
  if(eventColumnsAsFactor) {
    dataSetLong <- dataSetLong[, -which(names(dataSetLong) %in% eventColumns)]
  }
  
  # Combine results and output
  dataSetLong <- cbind(obj, timeInt, responses, dataSetLong)
  return(dataSetLong)  
}

##################################################################################
# Function for dataLong in the case of competing risks models 
# with time dependent covariates
##################################################################################

dataLongCompRisksTimeDep <- function (dataSet, timeColumn, eventColumns, 
                                      eventColumnsAsFactor=FALSE, idColumn, 
                                      timeAsFactor=TRUE) {
  
  # Alternative input for competing risks coding 
  # as factor instead of multiple binary columns
  if(eventColumnsAsFactor) {
    respFact <- dataSet[, eventColumns]
    respMat <- model.matrix(~., data.frame(respFact))[, -1]
    eventColumns <- levels(respFact)[-1]
    dimnames(respMat) [[2]] <- eventColumns
    dataSet <- cbind(respMat, dataSet)
  }
  
  # Construct censoring variable in short format
  indizeseventColumns <- sapply(1:length(eventColumns), function (x) which(eval(eventColumns [x])==colnames(dataSet)))
  NoeventColumns <- length(indizeseventColumns)
  # dataSet_eventColumns <- dataSet [, indizeseventColumns]
  # eventColumnsorShort <- 1 - rowSums(dataSet_eventColumns)
  responsesShort <- cbind(1, matrix(0, nrow=nrow(dataSet), ncol=NoeventColumns))
  namesEventColumnsNew <- paste("e", 0:NoeventColumns, sep="")
  dimnames(responsesShort) [[2]] <- namesEventColumnsNew
  dataSet <- cbind(responsesShort, dataSet)
  
  # Rearrange original data, that ID is in increasing order
  dataSet <- dataSet[order(dataSet [, idColumn]), ]
  
  # Split data by persons
  splitDataSet <- split(x=dataSet, f=dataSet [, idColumn])
  lengthSplitDataSet <- 1:length(splitDataSet)
  
  # Get count of replicates in each split
  Counts <- lapply(splitDataSet, function (x) c(diff(x [, timeColumn]), 1))
  SumCounts <- as.numeric(sapply(Counts, function (x) sum(x)))
  
  # Replicate indices in each split
  Indizes <- lapply(lengthSplitDataSet, function (x) rep(x=1:dim(splitDataSet [[x]])[1], times=Counts [[x]]))
  
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
    if(rowSumsPre[ cumSumSumCounts[i] ]==1) {
      # Replace censoring column with zero, if an event occurs
      dataList[cumSumSumCounts[i], namesEventColumnsNew[1] ] <- 0
      
      # Set respective event column to one
      relCol <- which(dataList [cumSumSumCounts[i], eventColumns]==1) + 1
      dataList[cumSumSumCounts[i], namesEventColumnsNew[relCol] ] <- 1
    }
  }
  
  # Remove unnecessary columns
  if(eventColumnsAsFactor) {
    dataList <- dataList[, -which(names(dataList) %in% eventColumns)]
  }
  
  # Combine results and output
  dataComplete <- cbind(obj=dataListID, timeInt=dataListTimeInt, dataList)
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
# timeColumn: Variable name of discrete time interval (character scalar)

# Output
# Original data.frame dataSetLong, but with added censoring process as first variable "yCens"

dataCensoring <- function(dataSetLong, respColumn, timeColumn){
  # Select last observed time
  lastTimeInd <- c((which(dataSetLong[, timeColumn]==1)-1)[-1], dim(dataSetLong)[1])
  
  # Replicate original y
  ResultVec <- dataSetLong[, respColumn]
  
  # Change relevant values of y to yCens
  ResultVec[lastTimeInd] <- ifelse(dataSetLong[lastTimeInd, "y"]==0, 1, NA)
  
  # Append results to output
  NewDataSet <- cbind(yCens = ResultVec, dataSetLong)
  return(NewDataSet)
}

# Slow, old version
# dataCensoring <- function (dataSetLong, respColumn, idColumn) {
# 
#   # Input checks
#   
#   if(!is.data.frame(dataSetLong)) {stop("Argument *dataSetLong* is not in the correct format! Please specify as data.frame object.")}
#   if(!(is.character(respColumn))) {stop("Argument *respColumn* is not in the correct format! Please specify as character.")}
#   if(length(respColumn)!=1) {stop("Argument *respColumn* is not in the correct format! Please specify as scalar with length one.")}
#   if(!(is.character(idColumn))) {stop("Argument *idColumn* is not in the correct format! Please specify as character.")}
#   if(length(idColumn)!=1) {stop("Argument *idColumn* is not in the correct format! Please specify as scalar with length one.")}
#   
#   # Can *respColumn* be accessed in *dataSetLong*?
#   if(any(class(tryCatch(dataSetLong [, respColumn], error=function (e) e)) == "error")) {
#     stop("*respColumn* is not available in *dataSetLong*! Please specify the correct column of observed times.")
#   }
#   # Can *idColumn* be accessed in *dataSetLong*?
#   if(any(class(tryCatch(dataSetLong [, idColumn], error=function (e) e)) == "error")) {
#     stop("*idColumn* is not available in *dataSetLong*! Please specify the correct column of the identification number.")
#   }
# 
#   # Check if response has only zero or one values
#   if(!all(as.numeric(as.character(dataSetLong [, respColumn]))==0 | as.numeric(as.character(dataSetLong [, respColumn]))==1)) {
#     stop("*idColumn* is not a binary vector! Please check, that events equals 1 and 0 otherwise.")
#   }
#   
#   ############
#   # Main
#   
#   # Splits data by persons
#   SplitData <- split(x=dataSetLong [, respColumn], f=dataSetLong [, idColumn])
#   
#   for(i in 1:length(SplitData)) {
#     
#     if(tail(SplitData [[i]], n=1)==0) {
#       SplitData [[i]] [length(SplitData [[i]])] <- 1
#     }
#     
#     else {
#       
#       if(tail(SplitData [[i]], n=1)==1) {
#         
#         if(length(SplitData [[i]])==1) {
#           SplitData [[i]] <- NA
#         }
#         
#         else {
#           SplitData [[i]] [length(SplitData [[i]])] <- NA
#         }
#         
#       }
#     }
#   }
#   
#   # Output
#   ResultVec <- do.call(c, SplitData)
#   NewDataSet <- cbind(yCens=ResultVec, dataSetLong)
#   return(NewDataSet)
# }

####################
# dataCensoringShort

dataCensoringShort <- function(dataSet, eventColumns, timeColumn){
  # Change relevant values of y to yCens
  ResultVec <- ifelse(rowSums(dataSet[, eventColumns, drop=FALSE])==0, 1, 0)

  # Exclude observations with time interval==1 and yCens==0
  deleteIndices <- ResultVec==0 & dataSet[, timeColumn]==1
  dataSet <- dataSet[!deleteIndices, ]
  ResultVec <- ResultVec[!deleteIndices]
  
  # Change time intervals for censoring process
  timeCens <- as.numeric(as.character(dataSet[, timeColumn]))
  timeCens[ResultVec==0] <- timeCens[ResultVec==0]-1

  # Append results to output
  NewDataSet <- cbind(yCens = ResultVec, timeCens=timeCens, dataSet)
  return(NewDataSet)
}

#################
# dataLongSubdist

dataLongSubDist <- function(dataSet, timeColumn, eventColumns, 
                             eventFocus, timeAsFactor=TRUE) {
  #######################
  # Construct long format
  
  # Reorder event columns
  # First column correspond to the event of focus
  eventColumns <- c(eventColumns[eventColumns==eventFocus], 
                    eventColumns[eventColumns!=eventFocus])
  
  # Construct object counter
  dataSet_timeColumn <- as.numeric(as.character(dataSet[, timeColumn]))
  tmax <- max(dataSet_timeColumn)
  delta_i <- rowSums(dataSet[, eventColumns])==1
  epsilon_i <- sapply(1:dim(dataSet)[1], 
                      function(x) ifelse(delta_i[x], 
                                         which(dataSet[x, eventColumns]==1), 0))
  dataSet_timeColumn <- ifelse(delta_i & 
                                 dataSet[, eventFocus] == 0, tmax, 
                               dataSet_timeColumn)
  obj <- rep(1:nrow(dataSet), each=tmax)
  
  # Construct time intervals
  timeInt <- c(sapply(1:length(dataSet_timeColumn), 
                      function(k) 1:tmax))
  
  # Construct response for event of interest
  eventFocus <- dataSet[, eventFocus]
  y <- c(sapply(1:length(dataSet_timeColumn), 
                function(k) c(rep(0, dataSet_timeColumn[k]-1), 
                              eventFocus[k], 
                              rep(0, tmax-dataSet_timeColumn[k]) ) ))
  
  # Estimation of weights
  estG <- estSurvCens(dataSet=dataSet, timeColumn=timeColumn, 
                      eventColumns=eventColumns)
  weights <- estG[timeInt] / 
    estG[pmin(dataSet[obj, timeColumn], timeInt)] *
    (ifelse(timeInt <= dataSet[obj, timeColumn], 1, 0) + 
       ifelse(dataSet[obj, timeColumn] <= (timeInt-1) &
                delta_i[obj] * epsilon_i[obj] > 1, 1, 0) )
  
  # Combine results
  if(timeAsFactor){timeInt <- factor(timeInt)}
  dataSetLong <- cbind(obj=obj, timeInt=timeInt, y=y, dataSet[obj, ])
  Output <- cbind(dataSetLong, subDistWeights=weights)
  return(Output)
} 

#################################
# Multi Spell data transformation

dataLongMultiSpell <- function(dataSet, timeColumn, censColumn, idColumn,
                               timeAsFactor=FALSE){
  
  # Expand responses
  # Assumption: Reference level corresponds to censoring
  if(!is.factor(dataSet[, censColumn])){
    dataSet[, censColumn] <- factor(dataSet[, censColumn])
  }
  eventExpand <- model.matrix(~.-1, dataSet[, censColumn, drop=FALSE])
  dimnames(eventExpand)[[2]] <- paste("e", 0:(dim(eventExpand)[2]-1), sep="")
  dataSet <- cbind(obj=NA, timeInt=NA, eventExpand, dataSet)
  
  # Extract indices for IDs from time
  splitDat <- split(dataSet, dataSet[, idColumn])
  lenSplitDat <- length(splitDat)
  
  # Find discrete time intervals with time gaps between spells
  for( j in 1:lenSplitDat ){
    
    # Calculate time differences
    actualTime <- splitDat[[j]][, timeColumn]
    diffTime <- diff( actualTime )
    # Define standard indices
    defIndices <- 1:length(actualTime)
    if(actualTime[1] > 1){
      defIndices <- c( defIndices, rep(1, actualTime[1]-1) )
    }
    
    for( i in 1:length(diffTime) ){
        if(length(diffTime) > 0){
          if(diffTime[i] > 1){
            defIndices <- c( defIndices, rep(i+1, diffTime[i]-1) )
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
      
      if(length(splitDat[[j]][defIndices==k, "e0"]) > 1){
        relInd <- which(defIndices==k)
        relInd <- relInd[-length(relInd)]
        splitDat[[j]][relInd, "e0"] <- 1
        splitDat[[j]][relInd, paste("e", 
                                    1:(dim(eventExpand)[2]-1), sep="")] <- 0
      }
      
    }
    
  }
  
  # Return full data set
  allDat <- as.data.frame(rbindlist(splitDat))
  
  # Should time be formatted as factor?
  if(timeAsFactor){
    allDat$timeInt <- factor(allDat$timeInt)
  }
  
  return(allDat)
}
