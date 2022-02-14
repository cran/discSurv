\name{evalIntPredErr}
\alias{evalIntPredErr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wrapper function to compute the integrated prediction error
}
\description{
Convenient version to directly compute the integrated prediction error curve without the need to use multiple functions.
}
\usage{
evalIntPredErr(hazPreds, survPreds=NULL, newTimeInput, newEventInput, 
trainTimeInput, trainEventInput, testDataLong, tmax=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hazPreds}{
  Predicted discrete hazards in the test data. The predictions have to be made with a data set in long format with prefix "dataLong", see e. g. \code{\link{dataLong}}, \code{\link{dataLongTimeDep}}.
}
\item{survPreds}{
Predicted discrete survival functions of all persons (list of numeric vectors). Each list element corresponds to the survival function of one person beginning with the first time interval. The default NULL assumes that arguments "hazPreds" and "testDataLong" are available. If the last two arguments are not available, then argument "survPreds" needs to be specified. 
}
\item{newTimeInput}{
Discrete time intervals in short format of the test set (integer vector).
}
\item{newEventInput}{
Events in short format in the test set (binary vector).
}
\item{trainTimeInput}{
Discrete time intervals in short format of the training set (integer vector).
}
\item{trainEventInput}{
Events in short format in the training set (binary vector).
}
\item{testDataLong}{
Test data in long format. Needs to be specified only, if argument "survPreds" is NULL. In this case the discrete survival function is calculated based on the predicted hazards. It is assumed that the data was preprocessed with a function with prefix "dataLong", see e. g. \code{\link{dataLong}}, \code{\link{dataLongTimeDep}}.
}
\item{tmax}{
Gives the maximum time interval for which prediction errors are calculated. It must be smaller than the maximum observed time in the training data of the object produced by function \code{\link{predErrDiscShort}}! Default=NULL means, that all observed intervals are used.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Integrated prediction error (numeric scalar).
}
\references{
Gerhard Tutz and Matthias Schmid, (2016), \emph{Modeling discrete time-to-event data}, Springer series in statistics, Doi: 10.1007/978-3-319-28158-2

Tilmann Gneiting and Adrian E. Raftery, (2007),
\emph{Strictly proper scoring rules, prediction, and
estimation},
Journal of the American Statistical Association 102 (477), 359-376
%@article {PredErrCurv,
%author={Tilmann Gneiting and Adrian E. Raftery},
%title={Strictly proper scoring rules, prediction, and
%estimation}, 
%journal={Journal of the American Statistical Association 102 (477), 359-376}, 
%year={2007}
%}
}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{predErrDiscShort}}, \code{\link{intPredErrDisc}}, \code{\link[stats]{aggregate}}
}
\examples{
##########################
# Example with cancer data

library(survival)
head(cancer)

# Data preparation and convertion to 30 intervals
cancerPrep <- cancer
cancerPrep$status <- cancerPrep$status-1
intLim <- quantile(cancerPrep$time, prob=seq(0, 1, length.out=30))
intLim [length(intLim)] <- intLim [length(intLim)] + 1

# Cut discrete time in smaller number of intervals
cancerPrep <- contToDisc(dataSet=cancerPrep, timeColumn="time", intervalLimits=intLim)

# Generate training and test data
set.seed(753)
TrainIndices <- sample (x=1:dim(cancerPrep) [1], size=dim(cancerPrep) [1]*0.75)
TrainCancer <- cancerPrep [TrainIndices, ]
TestCancer <- cancerPrep [-TrainIndices, ]
TrainCancer$timeDisc <- as.numeric(as.character(TrainCancer$timeDisc))
TestCancer$timeDisc <- as.numeric(as.character(TestCancer$timeDisc))

# Convert to long format
LongTrain <- dataLong(dataSet=TrainCancer, timeColumn="timeDisc", censColumn="status")
LongTest <- dataLong(dataSet=TestCancer, timeColumn="timeDisc", censColumn="status")
# Convert factors
LongTrain$timeInt <- as.numeric(as.character(LongTrain$timeInt))
LongTest$timeInt <- as.numeric(as.character(LongTest$timeInt))
LongTrain$sex <- factor(LongTrain$sex)
LongTest$sex <- factor(LongTest$sex)

# Estimate, for example, a generalized, additive model in discrete survival analysis
library(mgcv)
gamFit <- gam (formula=y ~ s(timeInt) + s(age) + sex + ph.ecog, data=LongTrain, family=binomial())
summary(gamFit)

# 1. Specification of predicted discrete hazards
# Estimate survival function of each person in the test data
testPredHaz <- predict(gamFit, newdata=LongTest, type="response")

# 1.1. Calculate integrated prediction error
evalIntPredErr(hazPreds=testPredHaz, survPreds=NULL, 
newTimeInput=TestCancer$timeDisc, newEventInput=TestCancer$status, 
trainTimeInput=TrainCancer$timeDisc, trainEventInput=TrainCancer$status, 
testDataLong=LongTest)

# 2. Alternative specification
oneMinusPredHaz <- 1-testPredHaz
predSurv <- aggregate(oneMinusPredHaz ~ obj, data=LongTest, FUN=cumprod, na.action=NULL)

# 2.1. Calculate integrated prediction error
evalIntPredErr(survPreds=predSurv[[2]], 
newTimeInput=TestCancer$timeDisc, newEventInput=TestCancer$status, 
trainTimeInput=TrainCancer$timeDisc, trainEventInput=TrainCancer$status)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line