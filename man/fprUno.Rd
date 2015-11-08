\name{fprUno}
\alias{fprUno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
False Positive Rate Uno
}
\description{
Estimates the predictive false positive rate (fpr) based on cross validation and generalized, linear models (see \code{\link{glm}}). The concept was suggested by Uno et al. (see references)
}
\usage{
fprUno(timepoint, dataSet, trainIndices, survModelFormula, 
censModelFormula, linkFunc = "logit", idColumn = NULL)
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
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
\itemize{
	\item{Output: }{List with objects:}
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
\code{\link{tprUno}}, \code{\link{tprUnoShort}}, \code{\link{aucUno}}, \code{\link{concorIndex}}, \code{\link[caret]{createFolds}}, \code{\link{glm}}
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

# Estimate false positive rate of time interval 7:
tryFPR <- fprUno (timepoint=7, dataSet=UnempDurSubset, trainIndices=CVfolds,  
survModelFormula=spell ~ age + logwage, censModelFormula=censor1 ~ 1, 
linkFunc="logit", idColumn=NULL)
tryFPR
plot(tryFPR)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line