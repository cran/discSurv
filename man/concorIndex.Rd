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
It is assumed that all time points up to the last interval [a_q, Inf) are available. If not already present, these can be added manually.
}

\seealso{
\code{\link{aucUno}}
}
\examples{

# Example with cross validation and unemployment data 
library(Ecdat)
library(caret)
data(UnempDur)

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

# Estimate false positive rate of time interval 7:
tryAUC <- aucUno (tprObj=tryTPR, fprObj=tryFPR)
tryAUC
plot(tryAUC)

\dontrun{# Estimate global concordance index:
tryConcorIndex <- concorIndex (tryAUC)
tryConcorIndex
summary(tryConcorIndex)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line