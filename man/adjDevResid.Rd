\name{adjDevResid}
\alias{adjDevResid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adjusted Deviance Residuals
}
\description{
Calculates the adjusted deviance residuals. This function only supports generalized, linear models. The adjusted deviance residuals should be approximately normal distributed, in the case of a well fitting model.
}
\usage{
adjDevResid(dataSet, survModelFormula, censColumn, linkFunc = "logit", idColumn = NULL)
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
		\item{AdjDevResid: } {Adjusted deviance residuals as numeric vector.}
		\item{GlmFit: } {Fit object of class (generalized, linear model used in the calculations)}
	}
	\item{Input: }{A list of given argument input values (saved for reference)}
}
}

\references{
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
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{devResid}}, \code{\link{brierScore}}, \code{\link{glm}}, \code{\link{predErrDiscShort}}
}
\examples{
# Example with cross validation and unemployment data 
library(Ecdat)
data(UnempDur)
summary(UnempDur$spell)

# Extract subset of data
set.seed(635)
IDsample <- sample(1:dim(UnempDur)[1], 100)
UnempDurSubset <- UnempDur [IDsample, ]

# Calculate adjusted deviance residuals for the unemployment data subset
adjDevianceResiduals <- adjDevResid (dataSet=UnempDurSubset, 
survModelFormula=spell ~ age + logwage, 
censColumn="censor1", linkFunc="logit", idColumn=NULL)
adjDevianceResiduals

# Exclude outliers
adjDevResidWoOut <- adjDevianceResiduals$Output$AdjDevResid [
adjDevianceResiduals$Output$AdjDevResid < 
quantile(adjDevianceResiduals$Output$AdjDevResid, prob=0.9) & 
adjDevianceResiduals$Output$AdjDevResid > 
quantile(adjDevianceResiduals$Output$AdjDevResid, prob=0.1)]

# Compare nonparametric density estimate of adjusted deviance residuals 
# with adapted normal distribution
plot(density(adjDevResidWoOut), xlim=c(-10,10), 
main="Density comparison: Normal vs nonparametric estimate", lwd=2, las=1, col="red")
dnorm1 <- function (x) {dnorm (x, mean=mean(adjDevResidWoOut), sd=sd(adjDevResidWoOut))}
curve (dnorm1, n=500, add=TRUE, col="black", lwd=2)
legend("topright", legend=c("Normal", "Nonpar"), lty=1, lwd=2, col=c("black", "red"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }