\name{brierScore}
\alias{brierScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Brier Score Calculation
}
\description{
Computes the brier score of each person based on a generalized, linear model with cross validation.
}
\usage{
brierScore(dataSet, trainIndices, survModelFormula, 
linkFunc = "logit", censColumn, idColumn = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
Original data in short format. Must be of class "data.frame".
}
  \item{trainIndices}{
Must be a list of indices used for training data sets. E. g. if a 10-fold cross validation should be conducted, then the list should have length 10 and in each element a integer vector of training indices.
}
  \item{survModelFormula}{
Gives the specified relationship of discrete response and covariates. The formula is designed, that the intercepts for the time dependent base line hazards are always included. Therefore only covariates should be given in this formula. This argument is required to be of class "formula".
}
  \item{linkFunc}{
Specifies the desired link function in use of generalized, linear models.
}
  \item{censColumn}{
Gives the column name of the event indicator (1=observed, 0=censored). Must be of type "character".
}
  \item{idColumn}{
Gives the column name of the identification number of each person. The argument must be of type "character". Default NULL means, that each row equals one person (no repeated measurements).
}
}
\details{
At normal circumstances a covariate free model for the censoring weights is sufficient.
}
\value{
Numeric vector of brier scores for each person in the original data.
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}

Matthias Schmid \email{matthias.schmid@imbie.uni-bonn.de}
}

\note{
It is assumed that all time points up to the last interval [a_q, Inf) are available in short format. If not already present, these can be added manually.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{adjDevResid}}, \code{\link{devResid}}, \code{\link{glm}}, \code{\link{predErrDiscShort}}
}
\examples{
# Example with cross validation and unemployment data 
library(Ecdat)
library(caret)
data(UnempDur)
summary(UnempDur$spell)

# Extract subset of data
set.seed(635)
IDsample <- sample(1:dim(UnempDur)[1], 100)
UnempDurSubset <- UnempDur [IDsample, ]
head(UnempDurSubset)
range(UnempDurSubset$spell)
set.seed(7550)
CVfolds <- createFolds (y=UnempDurSubset$spell, returnTrain=TRUE, k=2)

# Calculate brier score
tryBrierScore <- brierScore (dataSet=UnempDurSubset, trainIndices=CVfolds, 
survModelFormula=spell ~ age + logwage, linkFunc="logit", 
censColumn="censor1", idColumn=NULL)
tryBrierScore
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line