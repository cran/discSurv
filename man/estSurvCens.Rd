\name{estSurvCens}
\alias{estSurvCens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimated Survival Function
}
\description{
Estimates the marginal survival function G(T=t) of the censoring process. Compatible with single event and competing risks data.
}
\usage{
estSurvCens(dataSet, timeColumn, eventColumns)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
Data in original short format (data.frame).
}
  \item{timeColumn}{
Name of column with discrete time intervals (character scalar).
}
  \item{eventColumns}{
Names of the event columns of \code{dataSet}. In the competing risks case the event columns have to be in dummy encoding format (numeric vectors).  
}
}
%\details{
%
%}
\value{
Named vector of estimated survival function of the censoring process for all time points except the last theoretical interval. 
}
\references{
Gerhard Tutz and Matthias Schmid, (2016), \emph{Modeling discrete time-to-event data}, Springer series in statistics, Doi: 10.1007/978-3-319-28158-2
}

\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}

\note{
In the censoring survival function the last time interval [a_q, Inf) has the probability of zero.
}

\seealso{
\code{\link{estSurv}}
}
\examples{

# Load unemployment data
library(Ecdat)
data(UnempDur)

# Select subsample
subUnempDur <- UnempDur [1:100, ]

######################
# Single event example

# Estimate censoring survival function G(t)
estG <- estSurvCens(dataSet=subUnempDur, timeColumn="spell", 
eventColumns="censor1")
estG

#########################
# Competing risks example

# Estimate censoring survival function G(t)
estG <- estSurvCens(dataSet=subUnempDur, timeColumn="spell", 
eventColumns=c("censor1", "censor2", "censor3", "censor4"))
estG

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line