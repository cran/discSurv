\name{contToDisc}
\alias{contToDisc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Continuous to Discrete Transformation
}
\description{
Discretizes continuous time variable into a specified grid of censored data for discrete survival analysis. It is a data preprocessing step, before the data can be extendend in long format and further analysed with discrete survival models.  
}
\usage{
contToDisc(dataSet, timeColumn, intervalLimits, equi=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
Original data in short format. Must be of class "data.frame".
}
  \item{timeColumn}{
Name of the column with discrete survival times. Must be a scalar character value.
}
  \item{intervalLimits}{
Numeric vector of the right interval borders, e. g. if the intervals are [0, a_1), [a_1, a_2), [a_2, a_{max}), then intervalLimits = c(a_1, a_2, a_{max})
}
  \item{equi}{
Specifies if argument *intervalLimits* should be interpreted as number of equidistant intervals. Logical only TRUE or FALSE is allowed.
}
}
%\details{
%%
%}
\value{
Gives the data set expanded with a first column "timeDisc". This column includes the discrete time intervals (ordered factor). 
}
%\references{
%
%}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}
\note{
In discrete survival analysis the survival times have to be categorized in time intervals. Therefore this function is required, if there are observed continuous survival times.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dataLong}}, \code{\link{dataLongTimeDep}}, \code{\link{dataLongCompRisks}}
}
\examples{
# Example copenhagen stroke study data
library(pec)
data(cost)
head(cost)

# Convert observed times to months
# Right borders of intervals [0, a_1), [a_1, a_2), ... , [a_{\max-1}, a_{\max})
IntBorders <- 1:ceiling(max(cost$time)/30)*30

# Select subsample
subCost <- cost [1:100, ]
CostMonths <- contToDisc (dataSet=subCost, timeColumn="time", intervalLimits=IntBorders)
head(CostMonths)

# Select subsample giving number of equidistant intervals
CostMonths <- contToDisc (dataSet=subCost, timeColumn="time", intervalLimits=10, equi=TRUE)
head(CostMonths)
}
% Add one or more standard keywords, see file 'KEYWORDS' in % the R documentation directory.
\keyword{ datagen }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
