\name{predErrDiscShort}
\alias{predErrDiscShort}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction Error Curves for arbitrary prediction models
}
\description{
Estimates prediction error curves of arbitrary prediction models. In prediction error curves the estimated and observed survival functions are compared adjusted by weights at given timepoints.
}
\usage{
predErrDiscShort(timepoints, estSurvList, newTime, newEvent, trainTime, trainEvent)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{timepoints}{
Vector of the number of discrete time intervals. Must be of type integer.
}
  \item{estSurvList}{
List of persons in the test data. Each element contains a numeric vector of estimated survival functions of all given time points. 
}
  \item{newTime}{
Numeric vector of discrete survival times in the test data.
}
  \item{newEvent}{
Integer vector of univariate event indicator in the test data.
}
  \item{trainTime}{
Numeric vector of discrete survival times in the training data.
}
  \item{trainEvent}{
Integer vector of univariate event indicator in the training data.
}
}
\details{
The prediction error curves should be smaller than 0.25 for all time points, because this is equivalent to a random assignment error. 
}
\value{
\itemize{
	\item{List: }{List with objects:}
		\itemize{
			\item{Output: } {List with two components}
			\itemize{
				\item{predErr: } {Numeric vector with estimated prediction error values.
				Names give the evaluation time point.}
				\item{weights: } {List of weights used in the estimation. Each list component gives the weights of a person in the test data.}
			}
			\item{Input: } {A list of given argument input values (saved for reference)}
		}
	}
}
\references{
Van der Laan M. J. and J. M. Robins, (2003),
\emph{Unified Methods for Censored Longitudinal Data and Causality},
Springer, New York

Gerds T. A. and M. Schumacher, (2006),
\emph{Consistent estimation of the expected Brier score in general survival models with right-censored event times},
Biometrical Journal 48(6), 1029-1040
%@book{UniLongi,
%title={Unified Methods for Censored Longitudinal Data and Causality},
%author={Van der Laan, M. J. and J. M. Robins},
%publisher={Springer, New York},
%year={2003}
%}
%@article{ConsEstSurv,
%title={Consistent estimation of the expected Brier score in general survival models with right-censored event times},
%author={Gerds, T. A. and M. Schumacher},
%journal={Biometrical Journal 48(6), 1029-1040},
%year={2006}
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
\code{\link{aucUno}}, \code{\link[mgcv]{gam}}
}
\examples{
# Example with cross validation and unemployment data 
library(Ecdat)
library(mgcv)
data(UnempDur)
summary(UnempDur$spell)

# Extract subset of data
set.seed(635)
IDsample <- sample(1:dim(UnempDur)[1], 100)
UnempDurSubset <- UnempDur [IDsample, ]
head(UnempDurSubset)
range(UnempDurSubset$spell)

# Generate training and test data
set.seed(7550)
TrainIndices <- sample (x=1:dim(UnempDurSubset) [1], size=75)
TrainUnempDur <- UnempDurSubset [TrainIndices, ]
TestUnempDur <- UnempDurSubset [-TrainIndices, ]

# Convert to long format
LongTrain <- dataLong(dataSet=TrainUnempDur, timeColumn="spell", censColumn="censor1")
LongTest <- dataLong(dataSet=TestUnempDur, timeColumn="spell", censColumn="censor1")
# Convert factor to numeric for smoothing
LongTrain$timeInt <- as.numeric(as.character(LongTrain$timeInt))
LongTest$timeInt <- as.numeric(as.character(LongTest$timeInt))

######################################################################
# Estimate a generalized, additive model in discrete survival analysis

gamFit <- gam (formula=y ~ s(timeInt) + age + logwage, data=LongTrain, family=binomial())

# Estimate survival function of each person in the test data
oneMinusPredHaz <- 1 - predict(gamFit, newdata=LongTest, type="response")
predSurv <- aggregate(formula=oneMinusPredHaz ~ obj, data=LongTest, FUN=cumprod)

# Prediction error in first interval
tryPredErrDisc1 <- predErrDiscShort (timepoints=1, 
estSurvList=predSurv [[2]], newTime=TestUnempDur$spell,
newEvent=TestUnempDur$censor1, trainTime=TrainUnempDur$spell,
 trainEvent=TrainUnempDur$censor1)
tryPredErrDisc1
summary(tryPredErrDisc1)

# Prediction error of the 2. to 10. interval
tryPredErrDisc2 <- predErrDiscShort (timepoints=2:10,
estSurvList=predSurv [[2]], newTime=TestUnempDur$spell,
newEvent=TestUnempDur$censor1, trainTime=TrainUnempDur$spell,
trainEvent=TrainUnempDur$censor1)
tryPredErrDisc2
summary(tryPredErrDisc2)

########################################
# Fit a random discrete survival forest

library(randomForest)
LongTrainRF <- LongTrain
LongTrainRF$y <- factor(LongTrainRF$y)
rfFit <- randomForest (formula=y ~ timeInt + age + logwage, data=LongTrainRF)

# Estimate survival function of each person in the test data
oneMinusPredHaz <- 1 - predict(rfFit, newdata=LongTest, type="prob") [, 2]
predSurv <- aggregate(formula=oneMinusPredHaz ~ obj, data=LongTest, FUN=cumprod)

# Prediction error in first interval
tryPredErrDisc1 <- predErrDiscShort (timepoints=1, 
estSurvList=predSurv [[2]], newTime=TestUnempDur$spell,
newEvent=TestUnempDur$censor1, trainTime=TrainUnempDur$spell,
 trainEvent=TrainUnempDur$censor1)
tryPredErrDisc1
summary(tryPredErrDisc1)

# Prediction error of the 2. to 10. interval
tryPredErrDisc2 <- predErrDiscShort (timepoints=2:10,
estSurvList=predSurv [[2]], newTime=TestUnempDur$spell,
newEvent=TestUnempDur$censor1, trainTime=TrainUnempDur$spell,
trainEvent=TrainUnempDur$censor1)
tryPredErrDisc2
summary(tryPredErrDisc2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line