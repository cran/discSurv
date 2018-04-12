\name{adjDevResidShort}
\alias{adjDevResidShort}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adjusted Deviance Residuals in short format
}
\description{
Calculates the adjusted deviance residuals for arbitrary prediction models. The adjusted deviance residuals should be approximately normal distributed, in the case of a well fitting model.
}
\usage{
adjDevResidShort(dataSet, hazards)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
Data set in long format. Must be of class "data.frame".
}
  \item{hazards}{
Estimated hazard rates of the data in long format. Hazard rates are probabilities and therefore restricted to the interval [0, 1]
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
\itemize{
\item{Output: }{List with objects:}
	\itemize{
		\item{AdjDevResid: } {Adjusted deviance residuals as numeric vector}
	}
	\item{Input: }{A list of given argument input values (saved for reference)}
}
}

\references{
Gerhard Tutz and Matthias Schmid, (2016), \emph{Modeling discrete time-to-event data}, Springer series in statistics, Doi: 10.1007/978-3-319-28158-2

Gerhard Tutz, (2012), \emph{Regression for Categorical Data}, Cambridge University Press
%%@book {RegCat,
%%author = {Gerhard Tutz},
%%title = {Regression for Categorical Data},
%%publisher = {Cambridge University Press},
%%year = {2012}
%%}
}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}

\note{
The argument \emph{dataSet} must have a response with column name "y". The correct format of the dataset can be augmented by using \code{\link{dataLong}}. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{devResidShort}}, \code{\link{predErrDiscShort}}
}
\examples{
library(survival)

# Transform data to long format
heart[, "stop"] <- ceiling(heart[, "stop"])
set.seed(0)
Indizes <- sample(unique(heart$id), 25)
randSample <- heart[unlist(sapply(1:length(Indizes), 
function(x) which(heart$id==Indizes[x]))),]
heartLong <- dataLongTimeDep(dataSet=randSample, 
timeColumn="stop", censColumn="event", idColumn="id", timeAsFactor=FALSE)

# Fit a generalized, additive model and predict hazard rates on data in long format
library(mgcv)
gamFit <- gam(y ~ timeInt + surgery + transplant + s(age), data=heartLong, family="binomial")
hazPreds <- predict(gamFit, type="response")

# Calculate adjusted deviance residuals
devResiduals <- adjDevResidShort (dataSet=heartLong, hazards=hazPreds)$Output$AdjDevResid
devResiduals
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }