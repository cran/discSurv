\name{devResidShort}
\alias{devResidShort}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Deviance Residuals
}
\description{
Computes the root of the deviance residuals for evaluation of performance in discrete survival analysis.
}
\usage{
devResidShort(dataSet, hazards)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
Original data in short format. Must be of class "data.frame"
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
		\item{DevResid: } {Square root of deviance residuals as numeric vector.}
	}
	\item{Input: }{A list of given argument input values (saved for reference)}
}
}

\references{
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
\code{\link{adjDevResidShort}}, \code{\link{predErrDiscShort}}
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

# Calculate the deviance residuals
devResiduals <- devResidShort (dataSet=heartLong, hazards=hazPreds)$Output$DevResid

# Compare with estimated normal distribution
plot(density(devResiduals), 
main="Empirical density vs estimated normal distribution", 
las=1, ylim=c(0, 0.5))
tempFunc <- function (x) dnorm(x, mean=mean(devResiduals), sd=sd(devResiduals))
curve(tempFunc, xlim=c(-10, 10), add=TRUE, col="red")
# The empirical density seems like a mixture distribution, 
# but is not too far off in with values greater than 3 and less than 1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line