\name{concorIndex}
\alias{concorIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Concordance Index
}
\description{
Calculates the concordance index for discrete survival models (independent measure of time). This is the probability that, for a pair of randomly chosen comparable samples, the sample with the higher risk prediction will experience an event before the other sample or belongs to a higher binary class.
}
\usage{
concorIndex(aucObj, printTimePoints=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aucObj}{
Object of class "discSurvAucUno". This object is created using the function \code{\link{aucUno}}
}
  \item{printTimePoints}{
    Should messages be printed for each calculation of a discrete time interval? (logical scalar) Default is FALSE.
  }
}
\details{
The algorithm extracts all necessary information of the auc object (e. g. marginal probabilities and survival functions). 
}
\value{
List with objects
\itemize{
	\item{Output: }{Concordance index (named numeric vector)}
	\item{Input: }{List with all input arguments (saved for reference)}
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

\note{
It is assumed that all time points up to the last observed interval [a_{q-1}, a_q) are available. If not already present, these can be added manually.
}

\seealso{
\code{\link{aucUno}}
}
\examples{

# Example with cross validation and unemployment data 
library(Ecdat)
library(caret)
data(UnempDur)

# Evaluation of short term prediction for re-employed at full-time job
# Last interval q=14
# -> Set all time points with spell > 13 to time interval 13 and censored
lastObsInterval <- 13
UnempDurSubset <- UnempDur
UnempDurSubset[UnempDurSubset$spell > lastObsInterval, "censor1"] <- 0
UnempDurSubset[UnempDurSubset$spell > lastObsInterval, "spell"] <- lastObsInterval
head(UnempDurSubset)
range(UnempDurSubset$spell)

# Select cross-validation samples
set.seed(7550)
CVfolds <- createFolds (y=UnempDurSubset$spell, returnTrain=TRUE, k=2)

# Continuation ratio model formula
contModForm <- spell ~ logwage + ui + logwage*ui + age

# Estimate true positive rate of time interval 6: 
# Correspondes to three and a half month duration (each interval is of length two weeks)
tryTPR <- tprUno (timepoint=6, dataSet=UnempDurSubset, trainIndices=CVfolds, 
survModelFormula=contModForm, censModelFormula=censor1 ~ 1, 
linkFunc="logit", idColumn=NULL, timeAsFactor=FALSE)
tryTPR
plot(tryTPR)

# Estimate false positive rate of time interval 6:
tryFPR <- fprUno (timepoint=6, dataSet=UnempDurSubset, trainIndices=CVfolds,  
survModelFormula=contModForm, censModelFormula=censor1 ~ 1, 
linkFunc="logit", idColumn=NULL, timeAsFactor=FALSE)
tryFPR
plot(tryFPR)

# Estimate AUC rate of time interval 6:
tryAUC <- aucUno (tprObj=tryTPR, fprObj=tryFPR)
tryAUC
plot(tryAUC)

\dontrun{# Estimate global concordance index:
tryConcorIndex <- concorIndex (tryAUC, printTimePoints=TRUE)
tryConcorIndex
summary(tryConcorIndex)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line