#' Minimal Node Size Pruning
#'
#' Computes optimal minimal node size of a discrete survival tree from a given vector 
#' of possible node sizes by cross-validation. Laplace-smoothing can be applied to the 
#' estimated hazards.
#' 
#' @param formula Model formula for tree fitting("class formula")
#' @param data Discrete survival data in short format for which a survival tree is
#' to be fitted ("class data.frame").
#' @param treetype Type of tree to be fitted ("character vector"). Possible values are "rpart" or "ranger". The default
#' is to fit an rpart tree; when "ranger" is chosen, a ranger forest with a single tree is fitted.
#' @param splitruleranger String specifying the splitting rule of the ranger tree("character vector"). 
#' Possible values are either "gini", "extratrees" or "hellinger". Default is "hellinger".
#' @param sizes Vector of different node sizes to try ("integer vector"). 
#' Values should be non-negative.
#' @param indexList List of data partitioning indices for cross-validation ("class list").
#' Each element represents the test indices of one fold ("integer vector").
#' @param timeColumn Character giving the column name of the observed times in
#' the \emph{data} argument ("character vector").
#' @param eventColumn Character giving the column name of the event indicator in
#' the \emph{data} argument ("character vector").
#' @param lambda Parameter for laplace-smoothing. A value of 0 corresponds to 
#' no laplace-smoothing ("numeric vector").
#' @param logOut Logical value ("logical vector"). If the argument is set to TRUE, 
#' then computation progress will be written to console.
#' @details Computes the out-of-sample log likelihood for all data partitionings
#' for each node size in \emph{sizes} and returns the node size for which the log 
#' likelihood was minimal. Also returns an rpart tree with the optimal minimal 
#' node size using the entire data set.
#' @return A list containing the two items
#' \itemize{
#'   \item Optimal minimal node size - Node size with lowest out-of-sample log-likelihood
#'   \item tree - a tree object with type corresponding to \emph{treetype} argument with the optimal minimal node size
#' }
#' @examples
#' library(pec)
#' library(caret)
#' data(cost)
#' # Take subsample and convert time to years
#' cost$time <- ceiling(cost$time / 365)
#' costSub <- cost[1:50, ]
#' # Specify column names for data augmentation
#' timeColumn <- "time"
#' eventColumn <- "status"
#' # Create data partition for cross validation
#' indexList <- createFolds(costSub$status * max(costSub$time) + costSub$time, k = 5)
#' # specify function arguments and perform node size pruning
#' formula <- y ~ timeInt + prevStroke + age + sex
#' sizes <- 1:10
#' optiTree <- minNodePruning(formula, costSub, treetype = "rpart", sizes = sizes, 
#' indexList = indexList, timeColumn =  timeColumn, eventColumn = eventColumn, 
#' lambda = 1, logOut = TRUE)
#' @export minNodePruning
minNodePruning <- function(formula, data, treetype = "rpart", splitruleranger = "hellinger", sizes, indexList, timeColumn, 
                          eventColumn, lambda = 1, logOut = FALSE)
{
  #inputchecks
  if (!treetype %in% c("rpart", "ranger"))
  {
    stop("treetype must be either \"rpart\" or \"ranger\".")
  }
  mean_total_llh <- rep(NA, length(sizes))
  for (iNode in 1:length(sizes))
  {
    total_llh <- rep(NA, length(indexList))
    for (iTrainIndex in 1:length(indexList))
    {
      dataTrain <- data[-indexList[[iTrainIndex]],]
      dataTest <- data[indexList[[iTrainIndex]],]
      dataTrainLong <- dataLong(dataTrain, timeColumn, eventColumn)
      dataTestLong <- dataLong(dataTest, timeColumn, eventColumn)
      if(treetype == "ranger")
      {
        tree <- ranger(formula, dataTrainLong, num.trees = 1, mtry = length(attr(terms(formula), "term.labels")),
                      classification = TRUE, splitrule = splitruleranger, replace = FALSE, 
                      sample.fraction = 1, min.node.size = sizes[iNode])
        test_hazards <- survTreeLaplaceHazardRanger(tree, dataTrainLong, dataTestLong, lambda)
      } else
      {
        tree <- rpart(formula, dataTrainLong, method = "class", minbucket = sizes[iNode])
        test_hazards <- survTreeLaplaceHazard(tree, dataTestLong, lambda)
      }
      lh <- test_hazards[dataTestLong$y*nrow(test_hazards) + c(1:nrow(test_hazards))]
      llh <- -log(lh)
      total_llh[iTrainIndex] <- sum(llh)
    }
    mean_total_llh[iNode] <- mean(total_llh)
    if(logOut)
    {
      cat('\r',iNode/length(sizes)*100, "% finished")
      flush.console()
    }
  }
  optimalNodeSize <- sizes[which.min(mean_total_llh)]
  attr(optimalNodeSize,"llh") <- data.frame(sizes, mean_total_llh)
  dataLong <- dataLong(data, timeColumn, eventColumn)
  if(treetype == "ranger")
  {
    optimalTree <- ranger(formula, dataLong, num.trees = 1, mtry = length(attr(terms(formula), "term.labels")),
                         classification = TRUE, splitrule = "hellinger", replace = FALSE, 
                         sample.fraction = 1, min.node.size = optimalNodeSize)
  } else
  {
    optimalTree <- rpart(formula, dataLong, method = "class", minbucket = optimalNodeSize)
  }
  return(list("Optimal minimal node size" = optimalNodeSize,
              "tree" = optimalTree))
}
