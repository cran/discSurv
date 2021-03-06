\name{tprUno}
\alias{tprUno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
True positive Rate Uno
}
\description{
Estimates the true positive rate based on Uno et al. to evaluate the predictive accuracy of discrete generalized, linear survival models by cross validation.
}
\usage{
tprUno(timepoint, dataSet, trainIndices, survModelFormula, 
censModelFormula, linkFunc = "logit", idColumn = NULL, timeAsFactor=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{timepoint}{
Discrete time interval given that the false positive rate is evaluated (integer scalar)
}
  \item{dataSet}{
Original data in short format. Should be of class "data.frame".
}
  \item{trainIndices}{
List of Indices from original data used for training (list of integer vectors). The length of the list is equal to the number of cross validation samples.
}
  \item{survModelFormula}{
Formula of the discrete survival model. It is used in a generalized, linear model.
}
  \item{censModelFormula}{
Formula of the censoring model. It is used in a generalized, linear model. Usually this is done without covariates.
}
  \item{linkFunc}{
Link function of the generalized, linear model.
}
  \item{idColumn}{
Name of the column with identification numbers of persons. Default NULL means, that each row equals one person (no repeated measurements).
}
  \item{timeAsFactor}{
Should the time intervals be coded as factor? Default is to use factor. If the argument is false, the column is coded as numeric. 
}
}
\details{
  The formula \code{survModelFormula} must be in a specific structure: The response on the left side of the formula is the time of the short data format. On the right side are the covariates without time, e. g. Time ~ X1 + X2 if there are only two covariates. The time will be added automatically.
  
  The next formula \code{survModelFormula} has similiar structure. The difference is the left side of the formula: This must be the 
}
\value{
List with objects
\itemize{
	\item{Output: } {Data frame with two columns: "cutoff" gives the different marker values and "tpr" the true positive rates}
	\item{Input: } {A list of given argument input values (saved for reference). In addition there is the list element \code{orderMarker}, which gives the indices of the marker values in increasing order.}
}
}
\references{
Matthias Schmid, Gerhard Tutz and Thomas Welchowski, (2017), \emph{Discrimination Measures for Discrete Time-to-Event Predictions}, Econometrics and Statistics, Elsevier, 
Doi: 10.1016/j.ecosta.2017.03.008

Hajime Uno and Tianxi Cai and Lu Tian and L. J. Wei, (2007), \emph{Evaluating Prediction Rules for t-Year Survivors With Censored Regression Models}, Journal of the American Statistical Association

Patrick J. Heagerty and Yingye Zheng, (2005), \emph{Survival Model Predictive Accuracy and ROC Curves}, Biometrics 61, 92-105

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
\code{\link{tprUnoShort}}, \code{\link{fprUno}}, \code{\link{fprUnoShort}}, \code{\link{aucUno}}, \code{\link{concorIndex}}, \code{\link[caret]{createDataPartition}}, \code{\link{glm}}
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

# Estimate true positive rate of time interval 7: 
# Correspondes to three and a half month duration (each interval is of length two weeks)
tryTPR <- tprUno (timepoint=7, dataSet=UnempDurSubset, 
trainIndices=CVfolds, survModelFormula=spell ~ age + logwage, 
censModelFormula=censor1 ~ 1, linkFunc="logit", idColumn=NULL)
tryTPR
plot(tryTPR)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line