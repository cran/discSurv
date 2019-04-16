\name{evalCindex}
\alias{evalCindex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wrapper for evaluation of discrete concordance index
}
\description{
Convenient version to directly compute the discrete concordance without the need to use multiple functions.
}
\usage{
evalCindex(marker, newTime, newEvent, trainTime, trainEvent)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{marker}{
Gives the predicted values of the linear predictor of a regression model (numeric vector). May also be on the response scale. 
}
  \item{newTime}{
New time intervals in the test data (integer vector).
}
  \item{newEvent}{
New event indicators in the test data (integer vector with 0 or 1).
}
  \item{trainTime}{
Time intervals in the training data (integer vector).
}
  \item{trainEvent}{
Event indicators in the training data (integer vector with 0 or 1).
}
}
%\details{
%}
\value{
Value of discrete concordance index between zero and one (numeric scalar).
}
\references{
Matthias Schmid, Gerhard Tutz and Thomas Welchowski, (2017), \emph{Discrimination Measures for Discrete Time-to-Event Predictions}, Econometrics and Statistics, Elsevier, 
Doi: 10.1016/j.ecosta.2017.03.008

Hajime Uno and Tianxi Cai and Lu Tian and L. J. Wei, (2007), 
\emph{Evaluating Prediction Rules for t-Year Survivors With Censored Regression Models},
Journal of the American Statistical Association

Patrick J. Heagerty and Yingye Zheng, (2005), 
\emph{Survival Model Predictive Accuracy and ROC Curves}, 
Biometrics 61, 92-105
%@article {PredAcROC,
%author={Patrick J. Heagerty and Yingye Zheng},
%title={Survival Model Predictive Accuracy and ROC Curves}, 
%journal={Biometrics 61, 92-105}, 
%year={2005}
%}
%@article {EvalPredRules,
%author={Hajime Uno and Tianxi Cai and Lu Tian and L. J. Wei},
%title={Evaluating Prediction Rules for t-Year Survivors
%With Censored Regression Models}, 
%journal={Journal of the American Statistical Association}, 
%year={2007}
%}
}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}

\note{
It is assumed that all time points up to the last observed interval [a_{q-1}, a_q) are available.
}

\seealso{
\code{\link{tprUnoShort}}, \code{\link{fprUnoShort}}, \code{\link{aucUno}}, \code{\link{concorIndex}}
}
\examples{
##################################################
# Example with unemployment data and prior fitting

library(Ecdat)
library(caret)
library(mgcv)
data(UnempDur)
summary(UnempDur$spell)
# Extract subset of data
set.seed(635)
IDsample <- sample(1:dim(UnempDur)[1], 100)
UnempDurSubset <- UnempDur [IDsample, ]
set.seed(-570)
TrainingSample <- sample(1:100, 75)
UnempDurSubsetTrain <- UnempDurSubset [TrainingSample, ]
UnempDurSubsetTest <- UnempDurSubset [-TrainingSample, ]

# Convert to long format
UnempDurSubsetTrainLong <- dataLong(dataSet=UnempDurSubsetTrain, 
timeColumn="spell", censColumn="censor1")

# Estimate gam with smooth baseline
gamFit <- gam(formula=y ~ s(I(as.numeric(as.character(timeInt)))) + 
s(age) + s(logwage), data=UnempDurSubsetTrainLong, family=binomial())
gamFitPreds <- predict(gamFit, newdata=cbind(UnempDurSubsetTest, 
timeInt=UnempDurSubsetTest$spell))

# Evaluate C-Index based on short data format
evalCindex(marker=gamFitPreds, 
newTime=UnempDurSubsetTest$spell, 
newEvent=UnempDurSubsetTest$censor1, 
trainTime=UnempDurSubsetTrain$spell, 
trainEvent=UnempDurSubsetTrain$censor1)

#####################################
# Example National Wilm's Tumor Study

library(survival)
head(nwtco)
summary(nwtco$rel)

# Select subset
set.seed(-375)
Indices <- sample(1:dim(nwtco)[1], 500)
nwtcoSub <- nwtco [Indices, ]

# Convert time range to 30 intervals
intLim <- quantile(nwtcoSub$edrel, prob=seq(0, 1, length.out=30))
intLim [length(intLim)] <- intLim [length(intLim)] + 1
nwtcoSubTemp <- contToDisc(dataSet=nwtcoSub, timeColumn="edrel", intervalLimits=intLim)
nwtcoSubTemp$instit <- factor(nwtcoSubTemp$instit)
nwtcoSubTemp$histol <- factor(nwtcoSubTemp$histol)
nwtcoSubTemp$stage <- factor(nwtcoSubTemp$stage)

# Split in training and test sample
set.seed(-570)
TrainingSample <- sample(1:dim(nwtcoSubTemp)[1], round(dim(nwtcoSubTemp)[1]*0.75))
nwtcoSubTempTrain <- nwtcoSubTemp [TrainingSample, ]
nwtcoSubTempTest <- nwtcoSubTemp [-TrainingSample, ]

# Convert to long format
nwtcoSubTempTrainLong <- dataLong(dataSet=nwtcoSubTempTrain, 
timeColumn="timeDisc", censColumn="rel")

# Estimate glm
inputFormula <- y ~ timeInt + histol + instit + stage
glmFit <- glm(formula=inputFormula, data=nwtcoSubTempTrainLong, family=binomial())
linPreds <- predict(glmFit, newdata=cbind(nwtcoSubTempTest, 
timeInt=nwtcoSubTempTest$timeDisc))

# Evaluate C-Index based on short data format
evalCindex(marker=linPreds, 
newTime=as.numeric(as.character(nwtcoSubTempTest$timeDisc)), 
newEvent=nwtcoSubTempTest$rel, 
trainTime=as.numeric(as.character(nwtcoSubTempTrain$timeDisc)), 
trainEvent=nwtcoSubTempTrain$rel) 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line