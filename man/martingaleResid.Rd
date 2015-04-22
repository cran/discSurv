\name{martingaleResid}
\alias{martingaleResid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Martingale Residuals
}
\description{
Estimates the martingale residuals of a generalized, linear model.
}
\usage{
martingaleResid(dataSet, survModelFormula, censColumn, 
linkFunc = "logit", idColumn = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
Original data in short format. Should be of class "data.frame".
}
  \item{survModelFormula}{
Formula of the discrete survival model. It is used in a generalized, linear model.
}
  \item{censColumn}{
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
		\item{MartingaleResid: } {Square root of deviance residuals as numeric vector.}
		\item{GlmFit: } {Fit object of class (generalized, linear model used in the calculations)}
	}
	\item{Input: }{A list of given argument input values (saved for reference)}
}
}
\references{
Terry M. Therneau and Patricia M. Grambsch and Thomas R. Fleming, (1990),
\emph{Martingale-Based Residuals for Survival Models},
Biometrica, Vol. 77, No. 1, 147-160
%@article {MartResid,
%author={Terry M. Therneau and Patricia M. Grambsch and Thomas R. Fleming},
%title={Martingale-Based Residuals for Survival Models}, 
%journal={Biometrica, Vol. 77, No. 1, 147-160}, 
%year={1990}
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
\code{\link{tprUno}}, \code{\link{tprUnoShort}}, \code{\link{aucUno}}, \code{\link{concorIndex}}, \code{\link{glm}}
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

# Calculate martingale residuals for the unemployment data subset
MartResid <- martingaleResid (dataSet=UnempDurSubset, 
survModelFormula=spell ~ age + logwage, censColumn="censor1", 
linkFunc="logit", idColumn=NULL)
MartResid
sum(MartResid$Output$MartingaleResid)

# Plot martingale residuals vs each covariate in the event interval
# Dotted line is a loess estimate
plot(MartResid)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line