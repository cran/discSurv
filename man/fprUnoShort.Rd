\name{fprUnoShort}
\alias{fprUnoShort}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
False Positive Rate for arbitrary predition models
}
\description{
Estimates the false positive rate (based on concept of Uno, et al.) for an arbitrary discrete survival prediction model on one test data set. 
}
\usage{
fprUnoShort(timepoint, marker, newTime)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{timepoint}{
Gives the discrete time interval of which the fpr is evaluated (numeric scalar).
}
  \item{marker}{
Gives the predicted values of the linear predictor of a regression model (numeric vector). May also be on the response scale. 
}
  \item{newTime}{
New time intervals in the test data (integer vector).
}
}
\details{
This function is useful, if other models than generalized, linear models (glm) should be used for prediction. In the case of glm better use the cross validation version \code{\link{fprUno}}.
}
\value{
\itemize{
	\item{Output: }{A list with objects:}
		\itemize{
			\item{Output: } {Data frame with two columns: "cutoff" gives the different marker values and "fpr" the false positive rates}
			\item{Input: } {A list of given argument input values (saved for reference)}
		}
	}
}
\references{
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

Matthias Schmid \email{matthias.schmid@imbie.uni-bonn.de}
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{tprUno}}, \code{\link{tprUnoShort}}, \code{\link{aucUno}}, \code{\link{concorIndex}}, \code{\link{createFolds}}, \code{\link{glm}}
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
gamFitPreds <- predict(gamFit, newdata=cbind(UnempDurSubsetTest, timeInt=UnempDurSubsetTest$spell))

# Estimate false positive rate
fprGamFit <- fprUnoShort (timepoint=1, marker=gamFitPreds, newTime=UnempDurSubsetTest$spell)
plot(fprGamFit)

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

# Estimate tpr given one training and one test sample at time interval 10
fprFit <- fprUnoShort (timepoint=10, marker=linPreds, 
newTime=nwtcoSubTempTest$timeDisc)
plot(fprFit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line