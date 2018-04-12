\name{devResid}
\alias{devResid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Deviance Residuals
}
\description{
Computes the root of the deviance residuals for evaluation of performance in discrete survival analysis. A generalized, linear model is used for prediction.
}
\usage{
devResid(dataSet, survModelFormula, censColumn, linkFunc = "logit", idColumn = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
Original data in short format. Must be of class "data.frame".
}
  \item{survModelFormula}{
Gives the specified relationship of discrete response and covariates. The formula is designed, that the intercepts for the time dependent base line hazards are always included. Therefore only covariates should be given in this formula. This argument is required to be of class "formula".
}
  \item{censColumn}{
Gives the column name of the event indicator (1=observed, 0=censored). Must be of type "character".
}
  \item{linkFunc}{
Specifies the desired link function in use of generalized, linear models.
}
  \item{idColumn}{
Gives the column name of the identification number of each person. The argument must be of type "character". Default NULL means, that each row equals one person (no repeated measurements).
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
\itemize{
\item{Output: }{List with objects:}
	\itemize{
		\item{DevResid: } {Square root of deviance residuals as numeric vector.}
		\item{GlmFit: } {Fit object of class (generalized, linear model used in the calculations)}
	}
	\item{Input: }{A list of given argument input values (saved for reference)}
}
}

\references{
Gerhard Tutz and Matthias Schmid, (2016), \emph{Modeling discrete time-to-event data}, Springer series in statistics, Doi: 10.1007/978-3-319-28158-2

Gerhard Tutz, (2012), \emph{Regression for Categorical Data}, Cambridge University Press
%@book {RegCat,
%author = {Gerhard Tutz},
%title = {Regression for Categorical Data},
%publisher = {Cambridge University Press},
%year = {2012}
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
\code{\link{adjDevResid}}, \code{\link{brierScore}}, \code{\link{glm}}, \code{\link{predErrDiscShort}}
}
\examples{
library(Ecdat)
# Example with cross validation and unemployment data 
data(UnempDur)
summary(UnempDur$spell)

# Extract subset of data
set.seed(635)
IDsample <- sample(1:dim(UnempDur)[1], 100)
UnempDurSubset <- UnempDur [IDsample, ]

# Calculate deviance residuals for the unemployment data subset
devianceResiduals <- devResid (dataSet=UnempDurSubset, survModelFormula=spell ~ age + logwage, 
censColumn="censor1", linkFunc="logit", idColumn=NULL)
devianceResiduals
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line