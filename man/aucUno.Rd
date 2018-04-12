\name{aucUno}
\alias{aucUno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Area under the Curve Estimation
}
\description{
Estimates the time dependent area under the curve given calculated true positive rate and false positive rate. Both objects ("tprObj", "fprObj") and must have identical input arguments, e. g. same relationship of discrete response and covariates and supplied data sources. The values should be above 0.5 for a well predicting model, because random guessing would get this score. 
}
\usage{
aucUno(tprObj, fprObj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tprObj}{
Object of class "discSurvTprUno". See function \code{\link{tprUno}}
}
  \item{fprObj}{
Object of class "discSurvFprUno". See function \code{\link{fprUno}}
}
}
\details{
The auc is estimated by numerical integration of the ROC curve.
}
\value{
\itemize{
\item{Output: }{A list with objects:}
\itemize{
\item{Output: } {Named numeric vector with auc value}
\item{Input: } {A list of given argument input values (saved for reference)}
}
}
}

\references{
Matthias Schmid, Gerhard Tutz and Thomas Welchowski, (2017), \emph{Discrimination Measures for Discrete Time-to-Event Predictions, Econometrics and Statistics}, Elsevier, 
DOI: 10.1016/j.ecosta.2017.03.008

Hajime Uno and Tianxi Cai and Lu Tian and L. J. Wei, (2007), 
\emph{Evaluating Prediction Rules for t-Year Survivors With Censored Regression Models},
Journal of the American Statistical Association

Patrick J. Heagerty and Yingye Zheng, (2005), 
\emph{Survival Model Predictive Accuracy and ROC Curves}, 
Biometrics 61, 92-105
%%@article {PredAcROC,
%%author={Patrick J. Heagerty and Yingye Zheng},
%%title={Survival Model Predictive Accuracy and ROC Curves}, 
%%journal={Biometrics 61, 92-105}, 
%%year={2005}
%%}
%%@article {EvalPredRules,
%%author={Hajime Uno and Tianxi Cai and Lu Tian and L. J. Wei},
%%title={Evaluating Prediction Rules for t-Year Survivors
%%With Censored Regression Models}, 
%%journal={Journal of the American Statistical Association}, 
%%year={2007}
%}
}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}

Matthias Schmid \email{matthias.schmid@imbie.uni-bonn.de}
}
%%\note{
%%
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{tprUno}}, \code{\link{fprUno}}
}
\examples{
#####################################################
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

# Estimate true positive rate of time interval 7: 
# Correspondes to three and a half month duration (each interval is of length two weeks)
tryTPR <- tprUno (timepoint=7, dataSet=UnempDurSubset, trainIndices=CVfolds, 
survModelFormula=spell ~ age + logwage, censModelFormula=censor1 ~ 1, 
linkFunc="logit", idColumn=NULL, timeAsFactor=FALSE)
tryTPR
plot(tryTPR)

# Estimate false positive rate of time interval 7:
tryFPR <- fprUno (timepoint=7, dataSet=UnempDurSubset, trainIndices=CVfolds,  
survModelFormula=spell ~ age + logwage, censModelFormula=censor1 ~ 1, 
linkFunc="logit", idColumn=NULL, timeAsFactor=FALSE)
tryFPR
plot(tryFPR)

# Estimate auc
tryAUC <- aucUno (tprObj=tryTPR, fprObj=tryFPR)
tryAUC
plot(tryAUC)

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

# Estimate tpr given one training and one test sample at time interval 5
tprFit <- tprUnoShort (timepoint=5, marker=linPreds, 
newTime=nwtcoSubTempTest$timeDisc, newEvent=nwtcoSubTempTest$rel, 
trainTime=nwtcoSubTempTrain$timeDisc, trainEvent=nwtcoSubTempTrain$rel)

# Estimate fpr given one training and one test sample at time interval 5
fprFit <- fprUnoShort (timepoint=5, marker=linPreds, 
newTime=nwtcoSubTempTest$timeDisc)

# Estimate auc
tryAUC <- aucUno (tprObj=tprFit, fprObj=fprFit)
tryAUC
plot(tryAUC)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line